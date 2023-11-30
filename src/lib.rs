#![allow(dead_code)]

// Importer le module rand pour la génération de nombres aléatoires
use rand::*;
use std::f64::consts::PI;
use std::fmt;

// Point structure
#[derive(Default, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    fn new() -> Self {
        Default::default()
    }

    fn sqr_dist(self: Self, point: Point) -> f64 {
        (point.x - self.x).powi(2) + (point.y - self.y).powi(2)
    }
}

impl fmt::Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // HLSL style 'copy-paste'-able
        write!(f, "float2({}, {})", self.x, self.y)
    }
}

#[derive(Default)]
pub enum DomainType {
    #[default]
    Invalid,
    Box2D {
        w: f64,
        h: f64,
    },
    Disk {
        radius: f64,
    },
}

// Domain structure which hold a grid and the points on it
#[derive(Default)]
pub struct Domain {
    pub grid: Vec<Vec<usize>>,
    pub data: Vec<Point>,
    pub cell_size: f64,
    pub width: f64,
    pub height: f64,
    pub grid_width: usize,
    pub grid_height: usize,
    pub domain_type: DomainType,
}

impl Domain {
    pub fn box2d(width: f64, height: f64, sample_dist: f64) -> Self {
        const N: u8 = 2;

        let cell_size = sample_dist / (N as f64).sqrt();
        let ncells_width = (width / cell_size).ceil() as usize;
        let ncells_height = (height / cell_size).ceil() as usize;

        Domain {
            cell_size: cell_size,
            width: width,
            height: height,
            grid_width: ncells_width,
            grid_height: ncells_height,
            grid: vec![vec![usize::MAX; ncells_width]; ncells_height],
            domain_type: DomainType::Box2D { w: width, h: height },
            ..Default::default()
        }
    }

    pub fn disk(radius: f64, sample_dist: f64) -> Self {
        let mut domain = Domain {
            domain_type: DomainType::Disk { radius: radius },
            ..Default::default()
        };

        let n = domain.get_dimension();
        let cell_size = sample_dist / (n as f64).sqrt();
        let ncells = (2.0 * radius / cell_size).ceil() as usize;

        domain.cell_size = cell_size;
        domain.grid_width = ncells;
        domain.grid_height = ncells;
        domain.grid = vec![vec![usize::MAX; ncells]; ncells];

        domain
    }

    fn get_dimension(self: &Self) -> u8 {
        #![allow(unused_variables)]
        match self.domain_type {
            DomainType::Box2D { w, h } => 2,
            DomainType::Disk { radius } => 2,
            DomainType::Invalid => 0,
        }
    }

    fn extent(self: &Self) -> Point {
        match self.domain_type {
            DomainType::Box2D { w, h } => Point { x: w, y: h },
            DomainType::Disk { radius } => Point {
                x: 2.0 * radius,
                y: 2.0 * radius,
            },
            DomainType::Invalid => Point::default(),
        }
    }

    fn is_within_boundaries(self: &Self, point: Point) -> bool {
        match self.domain_type {
            DomainType::Box2D { w, h } => point.x.abs() < w / 2.0 && point.y.abs() < h / 2.0,
            DomainType::Disk { radius } => (point.x.powi(2) + point.y.powi(2)) < radius.powi(2),
            DomainType::Invalid => false,
        }
    }

    pub fn gen_random_point_within_domain(self: &Self) -> Point {
        match self.domain_type {
            DomainType::Box2D { w, h } => Point {
                x: rand::random::<f64>() * w - w / 2.0,
                y: rand::random::<f64>() * h - h / 2.0,
            },
            DomainType::Disk { radius } => {
                let sin_cos_theta = (rand::random::<f64>() * 2.0 * PI).sin_cos();
                let r = rand::random::<f64>() * radius;
                Point {
                    x: r * sin_cos_theta.0,
                    y: r * sin_cos_theta.1,
                }
            }
            DomainType::Invalid => Point::default(),
        }
    }

    pub fn insert_point(self: &mut Self, point: Point) {
        if !self.is_within_boundaries(point) {
            return;
        }

        let extent = self.extent();
        let x_index = ((point.x + extent.x / 2.0) / self.cell_size).floor() as usize;
        let y_index = ((point.y + extent.y / 2.0) / self.cell_size).floor() as usize;
        //assert!(self.grid[y_index][x_index] == usize::MAX);
        self.grid[y_index][x_index] = self.data.len();
        self.data.push(point);
    }

    pub fn is_valid_point(self: &Self, point: Point, radius: f64) -> bool {
        if !self.is_within_boundaries(point) {
            return false;
        }
        let extent = self.extent();
        let x_index = ((point.x + extent.x / 2.0) / self.cell_size).floor() as usize;
        let y_index = ((point.y + extent.y / 2.0) / self.cell_size).floor() as usize;

        let i0 = 0.max(x_index as i64 - 1) as usize;
        let i1 = (self.grid_width - 1).min(x_index + 1);
        let j0 = 0.max(y_index as i64 - 1) as usize;
        let j1 = (self.grid_height - 1).min(y_index + 1);

        // Rejects the point if too close to the others
        for j in j0..=j1 {
            for i in i0..=i1 {
                if self.grid[j][i] != usize::MAX && point.sqr_dist(self.data[self.grid[j][i]]) < radius.powi(2) {
                    return false;
                }
            }
        }

        return true;
    }
}

fn poisson_sample_internal(sample_dist: f64, k: usize, domain: &mut Domain) -> Vec<Point> {
    // The final list of points
    let mut points = Vec::<Point>::new();
    // The current "active" list of points
    let mut active = Vec::<Point>::new();

    // Initial random point p0
    let p0 = domain.gen_random_point_within_domain();

    domain.insert_point(p0);
    points.push(p0);
    active.push(p0);

    let mut rng = rand::thread_rng();

    while active.len() > 0 {
        // Pick a point 'p' from our active list
        let idx = rng.gen_range(0..active.len());
        let p_active = active[idx];

        // Try up to 'k' times to find a point that satifies:
        // - is at a distance between r and 2r from p
        // - is at a distance > r from nearby points
        let mut i = k;
        while i > 0 {
            let theta = rng.gen_range(0.0..2.0 * PI);

            let r = rng.gen_range(sample_dist..2.0 * sample_dist);

            let p = Point {
                x: p_active.x + r * theta.cos(),
                y: p_active.y + r * theta.sin(),
            };

            // If we succeed in finding a point, add to grid and list
            if domain.is_valid_point(p, sample_dist) {
                points.push(p);
                domain.insert_point(p);
                active.push(p);
                break;
            }

            i -= 1;
        }

        // Otherwise, remove point 'p' from the active list
        if i == 0 {
            active.remove(idx);
        }
    }

    points
}

pub fn poisson_sample_box2d(sample_dist: f64, k: usize, width: f64, height: f64) -> Vec<Point> {
    let mut domain = Domain::box2d(width, height, sample_dist);
    poisson_sample_internal(sample_dist, k, &mut domain)
}

pub fn poisson_sample_disk(sample_dist: f64, k: usize, domain_radius: f64) -> Vec<Point> {
    let mut domain = Domain::disk(domain_radius, sample_dist);
    poisson_sample_internal(sample_dist, k, &mut domain)
}
