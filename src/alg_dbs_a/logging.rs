#[cfg(feature = "logging")]
use std::cell::RefCell;

#[cfg(feature = "logging")]
#[derive(Default)]
pub struct Logs {
    pub loops: u64,
    pub cut_path_branches: u64,
    pub cut_leaf_branches: u64,
    pub full_reductions: u64,
    pub total_loops: u64,
    pub max_cache_size: u64,
    pub n_cache_flushes: u64,
}

/// Stores counters used to output information about algorithm runs
#[cfg(feature = "logging")]
impl Logs {
    pub fn print_logs_after_k(&self, k: usize) {
        for (name, val) in [
            ("total_loops", self.total_loops),
            ("loops", self.loops),
            ("cut_path_branches", self.cut_path_branches),
            ("cut_leaf_branches", self.cut_leaf_branches),
            ("full_reductions", self.full_reductions),
            ("max_cache_size", self.max_cache_size),
            ("n_cache_flushes", self.n_cache_flushes),
        ] {
            println!("#s k{}_{} {}", k, name, val);
        }

        if self.max_cache_size > 0 {
            for (name, val) in [
                ("max_cache_size", self.max_cache_size),
                ("n_cache_flushes", self.n_cache_flushes),
            ] {
                println!("#s k{}_{} {}", k, name, val);
            }
        }
    }

    pub fn reset_logs_for_new_k(&mut self) {
        self.loops = 0;
        self.cut_path_branches = 0;
        self.cut_leaf_branches = 0;
        self.full_reductions = 0;
        self.max_cache_size = 0;
        self.n_cache_flushes = 0;
    }
}

#[cfg(feature = "logging")]
thread_local! {
    pub static LOGS: RefCell<Logs> = RefCell::new(Logs::default());
}

#[cfg(feature = "logging")]
macro_rules! log {
    ($body:expr) => {
        LOGS.with(|s| $body(&mut *s.borrow_mut()))
    };
}

#[cfg(not(feature = "logging"))]
macro_rules! log {
    ($body:expr) => {};
}

pub(crate) use log;
