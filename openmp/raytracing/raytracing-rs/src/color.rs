use std::ops::{AddAssign, Mul, MulAssign};

#[derive(Debug, Copy, Clone)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Color {
    pub fn new(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b }
    }

    pub fn elementwise(self, rhs: Color) -> Color {
        Color {
            r: ((self.r as u32 * rhs.r as u32) / 255) as u8,
            g: ((self.g as u32 * rhs.g as u32) / 255) as u8,
            b: ((self.b as u32 * rhs.b as u32) / 255) as u8,
        }
    }

    pub fn black() -> Color {
        Color::new(0, 0, 0)
    }

    pub fn white() -> Color {
        Color::new(255, 255, 255)
    }
}

impl Mul<f64> for Color {
    type Output = Color;

    fn mul(self, rhs: f64) -> Color {
        Color {
            r: (self.r as f64 * rhs) as u8,
            g: (self.g as f64 * rhs) as u8,
            b: (self.b as f64 * rhs) as u8,
        }
    }
}

impl MulAssign<f64> for Color {
    fn mul_assign(&mut self, rhs: f64) {
        self.r = (self.r as f64 * rhs) as u8;
        self.g = (self.g as f64 * rhs) as u8;
        self.b = (self.b as f64 * rhs) as u8;
    }
}

impl AddAssign<Color> for Color {
    fn add_assign(&mut self, rhs: Color) {
        self.r += rhs.r;
        self.g += rhs.g;
        self.b += rhs.b;
    }
}
