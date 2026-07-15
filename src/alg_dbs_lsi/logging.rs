#[cfg(feature = "logging")]
use std::cell::RefCell;

#[cfg(feature = "logging")]
#[derive(Default)]
pub struct Logs {
    pub maf_calls: u64,
    pub br_lsi_calls: u64,
    pub br_process_calls: u64,
    pub br2_process_calls: u64,
    pub br_2_1_applied: u64,
    pub rr_2_2_1_applied: u64,
    pub br_2_2_2_applied: u64,
    pub br_2_2_2_special_case: u64,
    pub early_returns_br2_process_after_r1: u64,
    pub simple_lb_pruned: u64,
    pub lsi_cluster_reduc: u64,
    pub subtree_cluster_reduc: u64,
    pub subinstances_created: u64,
    pub subinstances_solved: u64,
    pub total_calls: u64,
}

#[cfg(feature = "logging")]
impl Logs {
    pub fn print_logs_after_k(&self, k: usize) {
        for (name, val) in [
            ("maf_calls", self.maf_calls),
            ("br_lsi_calls", self.br_lsi_calls),
            ("br_process_calls", self.br_process_calls),
            ("br2_process_calls", self.br2_process_calls),
            ("br_2_1_applied", self.br_2_1_applied),
            ("rr_2_2_1_applied", self.rr_2_2_1_applied),
            ("br_2_2_2_applied", self.br_2_2_2_applied),
            ("br_2_2_2_special_case", self.br_2_2_2_special_case),
            (
                "early_returns_br2_process_after_r1",
                self.early_returns_br2_process_after_r1,
            ),
            ("simple_lb_pruned", self.simple_lb_pruned),
            ("lsi_cluster_reduction", self.lsi_cluster_reduc),
            ("subinstances_created", self.subinstances_created),
            ("subinstances_solved", self.subinstances_solved),
            ("subtree_cluster_reduction", self.subtree_cluster_reduc),
        ] {
            println!("#s k{}_{} {}", k, name, val);
        }
    }
    pub fn print_logs_full_run(&self) {
        println!("#s total_calls {}", self.total_calls);
    }
    pub fn reset_logs_for_new_k(&mut self) {
        let temp_total = self.total_calls;
        *self = Logs::default();
        self.total_calls = temp_total;
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
