//#![feature(vec_into_raw_parts)]
///contains diophantine approximation to a rational number
mod dioph_appr;
pub use dioph_appr::dioph_appr;

///contains small pieces of original algorithm used in parallelized versions
pub mod par_pieces;
pub use par_pieces::b_seq_piece;
pub use par_pieces::s_piece;

///contains several subalgorithm for sieving in an interval
mod simple_sieve;
pub use simple_sieve::seg_sieve;
pub use simple_sieve::simple_seg_sieve;
pub use simple_sieve::simple_sieve;
pub use simple_sieve::sub_seg_sieve;

///original code from helfgott based on loops
mod sieve_loop;
pub use sieve_loop::sieve_loop;

///original code from helfgott based on iterators
mod sieve_s_iter;
pub use sieve_s_iter::sieve_s_iter;

///contains parallelization on the smalles iterator
mod s_par_sieve;
pub use s_par_sieve::s_par_sieve;

///contains parallelization on the biggest iterator
pub mod par_sieve_atomic;
pub use par_sieve_atomic::par_sieve_atomic;

//contains atomic parallel par_pieces WRONG
mod par_pieces_atomic;
pub use par_pieces_atomic::b_seq_piece_atomic;
pub use par_pieces_atomic::sub_seg_sieve_atomic;

//ptr version of paralleisation //TODO
//mod par_sieve_ptr;
//pub use par_sieve_ptr::par_sieve_ptr;

//ptr sub file
mod par_pieces_ptr;
pub use par_pieces_ptr::par_ptr;
//discord version
mod par_sieve_disc;
pub use par_sieve_disc::par_sieve_disc;
pub use par_sieve_disc::parallel;

pub mod parallel;
pub use parallel::parallel_sieve;
