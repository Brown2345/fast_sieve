mod dioph_appr;
mod new_seg_sieve;
mod parallel_sieve;
mod simple_sieve;
mod small_sieve;

pub use dioph_appr::dioph_appr;

pub use simple_sieve::seg_sieve;
pub use simple_sieve::simple_seg_sieve;
pub use simple_sieve::simple_sieve;
pub use simple_sieve::sub_seg_sieve;

pub use new_seg_sieve::new_seg_sieve;

pub use parallel_sieve::parallel_sieve;

pub use small_sieve::small_sieve;
