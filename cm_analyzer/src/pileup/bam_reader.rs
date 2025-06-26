use crate::Result;
use rust_htslib::bam::{self, IndexedReader, Read};


pub fn read_from_bam(bam_path: &str, target_id: u32) -> Result<IndexedReader> {
    let mut bam = bam::IndexedReader::from_path(bam_path).expect("Failed to open BAM file");
    bam.set_threads(1)?;
    bam.fetch(target_id)?;
    Ok(bam)
}

pub fn bam_to_iterator<'a>(
    bam: &'a mut IndexedReader,
) -> Result<impl Iterator<Item = (u32, u32, u32)> + 'a> {
    let iter = bam.pileup().into_iter().map(|pileup| {
        let pileup = pileup.unwrap();
        (pileup.tid(), pileup.pos(), pileup.depth())
    });
    Ok(iter)
}
