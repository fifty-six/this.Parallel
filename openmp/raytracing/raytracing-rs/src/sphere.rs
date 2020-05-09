use crate::Color;
use crate::Vec3;

#[derive(Debug)]
pub struct Sphere {
    pub radius: f64,
    pub diffuseness: f64,
    pub center: Vec3,
    pub color: Color,
}

const EPSILON: f64 = 0.01;

impl Sphere {
    pub fn intersect(&self, origin: Vec3, dir: Vec3) -> Option<(&Self, Vec3, f64)> {
        let diff = origin - self.center;

        let b = 2. * dir.dot(diff);
        let c = diff.dot(diff) - (self.radius * self.radius);

        let mut discriminant = (b * b) - 4. * c;

        if discriminant < 0. {
            return None;
        }

        discriminant = discriminant.sqrt();

        let sub = (-b - discriminant) / 2.;
        let add = (-b + discriminant) / 2.;

        // No hits either way
        if sub < 0. && add < 0. {
            return None;
        }

        // If they're both > 0, then we take the min
        // Otherwise take the max to ignore the negative value
        let mut min = if sub > 0. && add > 0. {
            f64::min(sub, add)
        } else {
            f64::max(sub, add)
        };

        min -= EPSILON;

        Some((self, origin + dir * min, min))
    }

    pub fn lambert_factor(&self, lambert: f64) -> f64 {
        1.0 - ((1.0 - lambert) * self.diffuseness)
    }
}
