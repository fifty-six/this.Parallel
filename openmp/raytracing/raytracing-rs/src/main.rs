mod color;
mod sphere;
mod vec3;

use color::Color;
use sphere::Sphere;
use vec3::Vec3;

use std::fs::File;
use std::io::{BufWriter, Error, Write};
use std::path::Path;

const SHADOW: f64 = 0.7;
const REFLECT: f64 = 0.7;

fn write_file(grid: Vec<Vec<Color>>, f_name: &str) -> Result<(), Error> {
    let path = Path::new(f_name);

    let file = File::create(path)?;

    let mut writer = BufWriter::new(file);

    writer.write_all(b"P3\n")?;
    writer.write_all(format!("{} {}\n", grid.first().unwrap().len(), grid.len()).as_bytes())?;
    writer.write_all(b"255\n")?;

    for y_row in grid {
        for px in y_row {
            writer.write_all(format!("{} {} {}\n", px.r, px.g, px.b).as_bytes())?;
        }
    }

    Ok(())
}

fn objects() -> Vec<Sphere> {
    vec![
        // Light
        Sphere {
            radius: 0.2,
            center: Vec3 {
                x: 0.,
                y: 1.25,
                z: -0.5,
            },
            color: Color {
                r: 255,
                g: 255,
                b: 255,
            },
        },
        // Floor
        Sphere {
            radius: 20000.25,
            center: Vec3 {
                x: 0.5,
                y: -20000.,
                z: 0.5,
            },
            color: Color {
                r: 234,
                g: 21,
                b: 81,
            },
        },
        // Blue sphere in middle
        Sphere {
            radius: 0.25,
            center: Vec3 {
                x: 0.5,
                y: 0.5,
                z: 0.5,
            },
            color: Color { r: 0, g: 0, b: 255 },
        },
        // Green sphere on right
        Sphere {
            radius: 0.25,
            center: Vec3 {
                x: 1.0,
                y: 0.5,
                z: 1.0,
            },
            color: Color { r: 0, g: 255, b: 0 },
        },
        // Red sphere on left
        Sphere {
            radius: 0.5,
            center: Vec3 {
                x: 0.0,
                y: 0.75,
                z: 1.25,
            },
            color: Color { r: 255, g: 0, b: 0 },
        },
    ]
}

pub fn get_color(
    origin: Vec3,
    dir: Vec3,
    light: Vec3,
    spheres: &[Sphere],
    depth: i8,
    max_depth: i8,
) -> Option<Color> {
    let hit = spheres
        .iter()
        .filter_map(|x| x.intersect(origin, dir))
        .min_by(|(_, _, t1), (_, _, t2)| t1.partial_cmp(t2).unwrap());

    let (sphere, intersection, _) = hit?;

    if sphere.color.r == 255 && sphere.color.g == 255 && sphere.color.b == 255 {
        return Some(sphere.color);
    }

    let light_dir = intersection.ray_to(light);

    let mut c = sphere.color;

    let shadow = spheres
        .iter()
        .skip(1)
        .any(|x| x.intersect(intersection, light_dir).is_some());

    let gradient = (intersection - sphere.center).as_normal();

    // Reflection
    {
        let reflect_ray = dir - gradient * 2. * gradient.dot(dir);

        let r_c = get_color(
            intersection,
            reflect_ray,
            light,
            spheres,
            depth + 1,
            max_depth,
        )
        .unwrap_or_else(Color::black);

        c *= 1. - REFLECT;

        c += r_c * REFLECT;
    }

    let mul = if shadow {
        // Binary shadow
        1. - SHADOW
    } else {
        // Lambert's law
        let cos_theta = f64::max(0., gradient.dot(light_dir));

        (1. - SHADOW) + SHADOW * cos_theta
    };

    c *= mul;

    Some(c)
}

fn main() {
    let max_y = 1080;
    let max_x = 1920;

    let eye = Vec3 {
        x: 0.5,
        y: 0.5,
        z: -1.,
    };

    let light = Vec3 {
        x: 0.,
        y: 1.25,
        z: -0.5,
    };

    let spheres = objects();

    let mut grid = vec![vec![Color::black(); max_x]; max_y];

    for (y, col) in grid.iter_mut().enumerate() {
        for (x, c) in col.iter_mut().enumerate() {
            let px_scaled = -0.25 + (x as f64 + 0.5) / max_y as f64;
            let py_scaled = ((max_y - y) as f64 + 0.5) / max_y as f64;

            let ray_dir = eye.ray_to(Vec3::new(px_scaled, py_scaled, 0.));

            if let Some(new_c) = get_color(eye, ray_dir, light, &spheres, 0, 10) {
                *c = new_c;
            }
        }
    }

    write_file(grid, "out.ppm").expect("Unable to write to file out.ppm!");
}
