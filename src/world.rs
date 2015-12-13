use std::io::prelude::*;
use std::fs::File;
use std::f64;
use std::marker::PhantomData;
use vec_map::VecMap;
use rand::{thread_rng, Rng};
use piston_window::*;
use gfx_device_gl::Resources;
use gfx::{self, Plane, Frame, Gamma};
use gfx::traits::{ToSlice, FactoryExt, Output, Stream, Factory};
use gfx::device::command::CommandBuffer;
use gfx::shade::TextureParam;
use gfx_graphics::GfxGraphics;
use nalgebra::{Pnt2, Vec2, Orig, FloatPnt, zero};
use acacia::{Tree, DataQuery, Node, AssociatedData, Position};
use acacia::partition::Ncube;

pub const WORLD_SIZE: f64 = 3000.0;
pub const GRAV: f64 = 100.0;
pub const TIME_SCALE: f64 = 12.0;
pub const EJECT_MASS: f64 = 0.05;
pub const EJECT_VEL: f64 = 30.0;
pub const COLLISION_EJECT_MASS: f64 = 0.05;
pub const COLLISION_EJECT_VEL: f64 = 15.0;
pub const NUM_BODIES: usize = 500;
//pub const SHADOW_LENGTH: f64 = 2000.0;

gfx_parameters!( ShadowParams {
    u_pos@ pos: [f32; 2],
    u_radius@ radius: f32,
    u_window@ window: [f32; 2],
    u_length@ length: f32,
    u_penumbra@ penumbra: TextureParam<R>,
});

gfx_vertex!( ShadowVertex {
    a_tex@ tex: [f32; 2],
});

gfx_parameters!( PostprocParams {
    u_tex@ tex: TextureParam<R>,
});

gfx_vertex!( PostprocVertex {
    a_pos@ pos: [f32; 3],
});

type GraphicsOutput<'a, O> = GfxGraphics<'a, Resources, CommandBuffer<Resources>, O>;

#[derive(Clone)]
pub struct Body {
    id: usize,
    x: f64,
    y: f64,
    dx: f64,
    dy: f64,
    mass: f64,
}

pub struct World {
    bodies: VecMap<Body>,
    next_body: usize,
    shadow_tex: Texture<Resources>,
    o_shadow: Frame<Resources>,
    o_postproc: Frame<Resources>,
    postproc: gfx::batch::Full<PostprocParams<Resources>>,
    //shadow: gfx::batch::Full<ShadowParams<Resources>>,
}

impl Position for Body {
    type Point = Pnt2<f64>;

    #[inline]
    fn position(&self) -> Pnt2<f64> {
        Pnt2::new(self.x, self.y)
    }
}


impl Body {
    pub fn new(id: usize) -> Body {
        let mut rng = thread_rng();
        Body {
            id: id,
            x: rng.gen_range(-WORLD_SIZE, WORLD_SIZE),
            y: rng.gen_range(-WORLD_SIZE, WORLD_SIZE),
            dx: rng.gen_range(-10.0, 10.0),
            dy: rng.gen_range(-10.0, 10.0),
            mass: rng.gen_range(2.0, 7.0)
        }
    }

    #[inline]
    pub fn radius(&self) -> f64 {
        5.0 * self.mass.sqrt()
    }

    pub fn draw<O>(&self, c: Context,
                   g: &mut GraphicsOutput<O>,
                   draw_state: &DrawState,
                   larger: bool) where O: Output<Resources>
    {
        let radius = self.radius() - 1.0;
        let border_color = if larger {
            [0.9, 0.3, 0.2, 1.0]
        }
        else {
            [0.5, 0.4, 0.3, 1.0]
        };
        Ellipse::new([0.8, 0.7, 0.5, 1.0]).border(ellipse::Border {
                color: border_color,
                radius: 1.0
            }).draw([self.x - radius, self.y - radius,
                     radius * 2.0, radius * 2.0],
                    &draw_state, c.transform, g);
    }

    pub fn draw_shadow<O>(&self,
                          c: Context,
                          g: &mut GraphicsOutput<O>,
                          draw_state: &DrawState,
                          view_x: f64,
                          view_y: f64,
                          view_width: f64,
                          view_height: f64) where O: Output<Resources>
    {
        //culling
        if self.x > view_x + view_width || self.y > view_y + view_height {
            return;
        }
        if self.x < view_x || self.y < view_y {
            if ((view_x - self.x) - (view_y - self.y)).abs() > view_width {
                return;
            }
        }

        let num = (self.radius() / 3.0).ceil().min(20.0) as i32;
        for angle in -num..num {
            let angle = -f64::consts::PI / 4.0 + angle as f64 / (64.0 * num as f64);
            let sx0 = self.x + self.radius() * angle.cos();
            let sy0 = self.y + self.radius() * angle.sin();
            let sx1 = self.x - self.radius() * angle.cos();
            let sy1 = self.y - self.radius() * angle.sin();
            let len = WORLD_SIZE * 10.0;

            Polygon::new([0.98, 0.98, 0.98, 1.0]).draw(
                &[[sx0, sy0], [sx1, sy1],
                  [sx1 + len * angle.cos(), sy1 - len * angle.sin()],
                  [sx0 + len * angle.cos(), sy0 - len * angle.sin()]],
                &draw_state.blend(draw_state::BlendPreset::Multiply), c.transform, g);
        }
    }
}

impl World {
    pub fn new(window: &PistonWindow, width: u32, height: u32) -> World {
        let mut bodies = VecMap::new();
        for i in 0..NUM_BODIES {
            bodies.insert(i, Body::new(i));
        }
        bodies[0].mass = 10.0;
        bodies[0].x = 0.0;
        bodies[1].y = 0.0;

        let mut tex_mem = Vec::with_capacity((width * height * 4) as usize);
        tex_mem.resize((width * height * 4) as usize, 0);
        let factory = &mut *window.factory.borrow_mut();

        let shadow_tex = Texture::from_memory_alpha(
            factory, &tex_mem[..],
            width, height, &TextureSettings::new()).unwrap();
        let o_shadow = Frame {
            width: width as u16,
            height: height as u16,
            colors: vec![Plane::Texture(shadow_tex.handle(), 0, None)],
            depth: None,
            stencil: None,
            gamma: Gamma::Original
        };

        let postproc_tex = Texture::from_memory_alpha(
            factory, &tex_mem[..],
            width, height, &TextureSettings::new()).unwrap();
        let o_postproc = Frame {
            width: width as u16,
            height: height as u16,
            colors: vec![Plane::Texture(postproc_tex.handle(), 0, None)],
            depth: None,
            stencil: None,
            gamma: Gamma::Original
        };

        let postproc = {
            let vertex_data = [
                PostprocVertex { pos: [-1.0, -1.0, 0.0] },
                PostprocVertex { pos: [1.0, -1.0, 0.0] },
                PostprocVertex { pos: [1.0, 1.0, 0.0] },
                PostprocVertex { pos: [-1.0, -1.0, 0.0] },
                PostprocVertex { pos: [1.0, 1.0, 0.0] },
                PostprocVertex { pos: [-1.0, 1.0, 0.0] }
            ];
            let mesh = factory.create_mesh(&vertex_data);
            let slice = mesh.to_slice(gfx::PrimitiveType::TriangleList);
            let mut vertex_src = Vec::new();
            let mut fragment_src = Vec::new();
            File::open("shaders/postproc_v.glsl").unwrap()
                .read_to_end(&mut vertex_src).unwrap();
            File::open("shaders/postproc_f.glsl").unwrap()
                .read_to_end(&mut fragment_src).unwrap();
            let program = factory.link_program(&vertex_src[..], &fragment_src[..]).unwrap();
            let state = gfx::DrawState::new();
            let sampler = factory.create_sampler(
                gfx::tex::SamplerInfo::new(gfx::tex::FilterMethod::Scale,
                                           gfx::tex::WrapMode::Clamp));
            let data = PostprocParams {
                tex: (postproc_tex.handle(), Some(sampler)),
                _r: PhantomData
            };

            let mut batch = gfx::batch::Full::new(mesh, program, data).unwrap();
            batch.slice = slice;
            batch.state = state;
            batch
        };

        /*let shadow = {
            let vertex_data = [
                ShadowVertex { tex: [1.0, 0.0] },
                ShadowVertex { tex: [1.0, 1.0] },
                ShadowVertex { tex: [0.0, 1.0] },
                ShadowVertex { tex: [0.0, 1.0] },
                ShadowVertex { tex: [0.0, 0.0] },
                ShadowVertex { tex: [1.0, 0.0] },
            ];
            let mesh = factory.create_mesh(&vertex_data);
            let slice = mesh.to_slice(gfx::PrimitiveType::TriangleList);
            let mut vertex_src = Vec::new();
            let mut fragment_src = Vec::new();
            File::open("shaders/shadow_v.glsl").unwrap()
                .read_to_end(&mut vertex_src).unwrap();
            File::open("shaders/shadow_f.glsl").unwrap()
                .read_to_end(&mut fragment_src).unwrap();
            let program = factory.link_program(&vertex_src[..], &fragment_src[..]).unwrap();
            let state = gfx::DrawState::new()
                .blend(gfx::BlendPreset::Multiply);

            let penumbra_tex = Texture::from_path(
                factory, "assets/penumbra.png", Flip::None, &TextureSettings::new()).unwrap();

            let data = ShadowParams {
                pos: [0.0, 0.0],
                radius: 0.0,
                window: [width as f32, height as f32],
                length: 0.0,
                penumbra: (penumbra_tex.handle(), None),
                _r: PhantomData
            };

            let mut batch = gfx::batch::Full::new(mesh, program, data).unwrap();
            batch.slice = slice;
            batch.state = state;
            batch
        };*/

        World {
            bodies: bodies,
            next_body: NUM_BODIES,
            shadow_tex: shadow_tex,
            o_shadow: o_shadow,
            o_postproc: o_postproc,
            postproc: postproc,
            //shadow: shadow
        }
    }

    #[inline]
    pub fn space_scale(&self) -> f64 {
        20.0 / self.bodies[0].radius()
    }

    pub fn eject_mouse(&mut self, x: f64, y: f64, width: f64, height: f64) {
        let x = (x - width / 2.0) / self.space_scale();
        let y = (y - height / 2.0) / self.space_scale();
        let dist = (x * x + y * y).sqrt();
        let dx = x / dist;
        let dy = y / dist;
        self.eject(0, EJECT_VEL, dx, dy, EJECT_MASS);
    }

    pub fn eject(&mut self, id: usize, vel: f64, dx: f64, dy: f64, mass: f64) {
        let mut body = Body {
            id: self.next_body,
            x: 0.0,
            y: 0.0,
            dx: dx * vel,
            dy: dy * vel,
            mass: mass * self.bodies[id].mass
        };
        {
            let body2 = &mut self.bodies[id];
            let dv = ((vel * vel * body.mass) / (body2.mass - body.mass)).sqrt();
            body.dx += body2.dx;
            body.dy += body2.dy;
            body2.dx -= dx * dv;
            body2.dy -= dy * dv;
            body2.mass -= mass * body2.mass;
            body.x = body2.x + (body2.radius() + body.radius() + 1.0) * dx;
            body.y = body2.y + (body2.radius() + body.radius() + 1.0) * dy;
        }
        self.bodies.insert(self.next_body, body);
        self.next_body += 1;
    }

    pub fn tick(&mut self, dt: f64) -> bool {
        let mut rng = thread_rng();

        let dt = dt * TIME_SCALE;

        let origin = Orig::orig();
        let bodies = self.bodies.clone();
        let low_x = bodies.values().fold(f64::MAX, |v, body| body.x.min(v));
        let low_y = bodies.values().fold(f64::MAX, |v, body| body.y.min(v));
        let high_x = bodies.values().fold(-f64::MAX, |v, body| body.x.max(v));
        let high_y = bodies.values().fold(-f64::MAX, |v, body| body.y.max(v));
        //println!("{}", (-low_x).max(-low_y).max(high_x).max(high_y));
        let tree = Tree::new(
            bodies.values(),
            Ncube::new(origin, 2.0 * (-low_x).max(-low_y).max(high_x).max(high_y) + 10.0),
            (origin, 0.0, 0),
            &|body: &&Body| (Pnt2::new(body.x, body.y), body.mass, body.id),
            &|&(pos1, mass1, _), &(pos2, mass2, _)| if mass1 + mass2 > 0.0 {
                (origin + (pos1.to_vec() * mass1 + pos2.to_vec() * mass2) / (mass1 + mass2),
                 mass1 + mass2, 0)
            }
            else {
                (origin, 0.0, 0)
            });

        let mut collisions = Vec::new();
        for (i, body) in self.bodies.iter_mut() {
            let ref_pos = Pnt2::new(body.x, body.y);
            let accel: Vec2<f64> = tree.query_data(|node| {
                let &(ref pos, _, _) = node.data();
                let d = FloatPnt::dist(&ref_pos, pos);
                let d2 = FloatPnt::dist(&node.partition().center(), pos);
                d < node.partition().width() / 0.5 + d2
            }).map(|&(pos, mass, id)| if id != i {
                    let dist = pos.dist(&ref_pos);
                    if dist < body.radius() + bodies[id].radius() && id != 0 {
                        collisions.push((i, id));
                    }
                    (pos - ref_pos) * (GRAV * mass / (dist * dist)) / dist
                }
                else {
                    zero()
                }).fold(zero(), |a, b| a + b);
            body.dx += dt * accel.x;
            body.dy += dt * accel.y;
            body.x += dt * body.dx;
            body.y += dt * body.dy;
        }

        for (a, b) in collisions.into_iter() {
            if self.bodies.contains_key(&a) && self.bodies.contains_key(&b) {
                let (dx, dy, mass) = {
                    let ba = &self.bodies[a];
                    let bb = &self.bodies[b];
                    let mass = ba.mass + bb.mass;
                    let dx = (ba.mass * ba.dx + bb.mass * bb.dx) / mass;
                    let dy = (ba.mass * ba.dy + bb.mass * bb.dy) / mass;
                    (dx, dy, mass)
                };
                let rel_dx = self.bodies[a].dx - self.bodies[b].dx;
                let rel_dy = self.bodies[a].dy - self.bodies[b].dy;
                let rel_v = (rel_dx * rel_dx + rel_dy * rel_dy).sqrt();
                let mass_ratio = self.bodies[a].mass / self.bodies[b].mass;
                //find the larger body
                let (new, mass_ratio) = if self.bodies[a].mass > self.bodies[b].mass {
                    if b == 0 {
                        return true;
                    }
                    self.bodies.remove(&b);
                    (a, mass_ratio)
                }
                else {
                    if a == 0 {
                        return true;
                    }
                    self.bodies.remove(&a);
                    (b, 1.0 / mass_ratio)
                };
                self.bodies[new].dx = dx;
                self.bodies[new].dy = dy;
                self.bodies[new].mass = mass;
                //check relative velocity
                if rel_v > 15.0 && mass_ratio < 4.0 {
                    for _ in 0..8 {
                        let angle = rng.gen_range(0.0, 2.0 * f64::consts::PI);
                        self.eject(new, COLLISION_EJECT_VEL,
                                   angle.cos(), angle.sin(),
                                   COLLISION_EJECT_MASS * rng.gen_range(0.0, 1.0));
                    }
                }
            }
        }

        false
    }

    pub fn draw(&mut self, window: &PistonWindow, stream: &mut GfxStream, width: u32, height: u32) {
        let c_abs = Context::new_abs(width as f64, height as f64);
        let view_x = -(width as f64 / 2.0 - self.bodies[0].x * self.space_scale());
        let view_y = -(height as f64 / 2.0 - self.bodies[0].y * self.space_scale());
        let c = c_abs
             .trans(-view_x, -view_y)
             .scale(self.space_scale(), self.space_scale());
        let draw_state = c.draw_state;//.multi_sample();

        {
            let renderer = &mut stream.ren;
            let g2d = &mut window.g2d.borrow_mut();
            {
                /*let mut stream = (&mut stream.ren, &self.o_shadow);
                stream.clear(ClearData {
                    color: [1.0, 1.0, 1.0, 0.0],
                    depth: 0.0,
                    stencil: 0
                });

                for body in self.bodies.values() {
                    self.shadow.params.pos = [
                        ((body.x - self.bodies[0].x) * self.space_scale()) as f32,
                        ((body.y - self.bodies[0].y) * self.space_scale()) as f32
                    ];
                    self.shadow.params.radius = (body.radius() * self.space_scale()) as f32;
                    self.shadow.params.length = (SHADOW_LENGTH * self.space_scale()) as f32;
                    stream.draw(&self.shadow).unwrap();
                }*/

                let g = &mut GfxGraphics::new(renderer, &self.o_shadow, g2d);
                clear([1.0, 1.0, 1.0, 1.0], g);
                for body in self.bodies.values() {
                    body.draw_shadow(c, g, &draw_state,
                                     view_x / self.space_scale(),
                                     view_y / self.space_scale(),
                                     width as f64 / self.space_scale(),
                                     height as f64 / self.space_scale());
                }
            }
            {
                let g = &mut GfxGraphics::new(renderer, &self.o_postproc, g2d);

                clear([0.05, 0.05, 0.15, 1.0], g);

                for body in self.bodies.values() {
                    body.draw(c, g, &draw_state, body.mass > self.bodies[0].mass);
                }

                let c = c_abs.scale(1.0, -1.0).trans(0.0, -(height as f64));
                Image::new().draw(&self.shadow_tex,
                                  &draw_state.blend(draw_state::BlendPreset::Multiply),
                                  c.transform, g);
            }
        }

        //post processing
        stream.draw(&self.postproc).unwrap();
    }
}
