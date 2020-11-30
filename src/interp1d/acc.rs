pub struct InterpAccel {
    pub(crate) cache: usize,
    pub(crate) hit_count: usize,
    pub(crate) miss_count: usize,
}

impl Default for InterpAccel {
    fn default() -> Self {
        Self::new()
    }
}

impl InterpAccel {
    pub fn new() -> InterpAccel {
        InterpAccel {
            cache: 0,
            hit_count: 0,
            miss_count: 0,
        }
    }
    pub fn reset(&mut self) {
        self.cache = 0;
        self.miss_count = 0;
        self.hit_count = 0;
    }
}
