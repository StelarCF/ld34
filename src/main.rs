extern crate piston_window;
extern crate rand;
extern crate vec_map;
extern crate gfx_device_gl;
#[macro_use]
extern crate gfx;
extern crate gfx_graphics;
extern crate acacia;
extern crate nalgebra;
extern crate music;
extern crate glutin;

use piston_window::*;

mod world;
use world::World;

enum GameState {
    Playing,
    Dead
}

fn main() {
    let (width, height) = glutin::get_available_monitors().next().map(|window|
        window.get_dimensions()).unwrap_or((1024, 768));

    //TODO configurable resolution
    let window: PistonWindow = WindowSettings::new("Protoplanet", (width, height))
        .exit_on_esc(true)
        .build()
        .unwrap();

    music::start::<(), _>(|| {
        music::bind_file((), "./assets/reverie.mp3");
        music::play(&(), music::Repeat::Forever);

        let mut world = World::new(&window, width, height);

        let mut paused = false;
        let mut mouse_pos = [0.0, 0.0];

        let mut glyphs = Glyphs::new(
            "./assets/OptimusPrinceps.ttf",
            window.factory.borrow().clone()).unwrap();

        let mut state = GameState::Playing;

        for e in window {
            e.update(|args| {
                match state {
                    GameState::Playing => if !paused {
                        if world.tick(args.dt) {
                            state = GameState::Dead;
                        }
                    },
                    GameState::Dead => ()
                }
            });

            e.draw_3d(|stream| {
                match state {
                    GameState::Playing => world.draw(&e, stream, width, height),
                    GameState::Dead => ()
                }
            });

            e.draw_2d(|c, g| {
                match state {
                    GameState::Dead => {
                        clear([0.05, 0.05, 0.15, 1.0], g);
                        Text::new_color([0.5, 0.2, 0.2, 1.0], 100)
                            .draw("You were absorbed...", &mut glyphs,
                                  &c.draw_state, c.trans(75.0, 100.0).transform, g);
                        Text::new_color([0.4, 0.4, 0.4, 1.0], 60)
                            .draw(&format!("Mass: {:.01}", world.bodies[0].mass)[..],
                                  &mut glyphs, &c.draw_state, c.trans(100.0, 300.0).transform, g);
                        Text::new_color([0.4, 0.4, 0.4, 1.0], 50)
                            .draw("Press any key to restart", &mut glyphs,
                                  &c.draw_state, c.trans(100.0, 500.0).transform, g);
                    },
                    GameState::Playing => if paused {
                        Text::new_color([0.5, 0.2, 0.2, 0.7], 50)
                            .draw("Paused", &mut glyphs,
                                &c.draw_state, c.trans(0.0, height as f64).transform, g);
                    }
                }
            });

            e.mouse_cursor(|x, y| {
                mouse_pos = [x, y];
            });

            e.press(|k| {
                match state {
                    GameState::Playing => match k {
                        Button::Mouse(MouseButton::Left) => world.eject_mouse(
                            mouse_pos[0], mouse_pos[1], width as f64, height as f64),
                        Button::Keyboard(Key::Space) => paused = !paused,
                        _ => ()
                    },
                    GameState::Dead => {
                        //restart game
                        state = GameState::Playing;
                        world = World::new(&e, width, height);
                    },
                }
            });
        }
    });
}
