use std::{
    fs::File,
    io::{BufWriter, Error, Write},
    panic,
    path::Path,
};

use crate::Color;

pub struct Canvas {
    pub width: usize,
    pub height: usize,
    data: Box<[u8]>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            data: vec![0; width * height * std::mem::size_of::<Color>()].into_boxed_slice(),
        }
    }

    pub fn set(&mut self, x: usize, y: usize, color: Color) {
        if x > self.width || y > self.height {
            panic!(
                "Out of bounds! x = {x}, y = {y} on container of size {width}, {height}",
                x = x,
                y = y,
                width = self.width,
                height = self.height
            );
        }

        let offset = (self.width * y + x) * std::mem::size_of::<Color>();

        self.data[offset] = color.r;
        self.data[offset + 1] = color.g;
        self.data[offset + 2] = color.b;
    }

    pub fn write(&self, path: &Path) -> Result<(), Error> {
        let file = File::create(path)?;

        let mut writer = BufWriter::new(file);

        writer.write_all(b"P6\n")?;
        writer.write_all(format!("{} {}\n", self.width, self.height).as_bytes())?;
        writer.write_all(b"255\n")?;

        writer.write_all(&self.data)?;

        Ok(())
    }
}
