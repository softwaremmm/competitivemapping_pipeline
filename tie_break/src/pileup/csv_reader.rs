use crate::Result;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::iter::Peekable;
use std::path::Path;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Aln {
    target_id: u32,
    ref_start: u32,
    ref_end: u32, // inclusive end
}

#[derive(Debug)]
pub struct CsvPileupIterator<I>
where
    I: Iterator<Item = Aln>,
{
    input: Peekable<I>,
    active: Vec<Aln>,
    current_target_id: Option<u32>,
    pos: u32,
}

impl<I> CsvPileupIterator<I>
where
    I: Iterator<Item = Aln>,
{
    pub fn new(input: I) -> Self {
        Self {
            input: input.peekable(),
            active: Vec::new(),
            current_target_id: None,
            pos: 0,
        }
    }
}

impl CsvPileupIterator<std::vec::IntoIter<Aln>> {
    pub fn from_path<P: AsRef<Path>>(path: P, target_id: Option<u32>) -> Result<Self> {
        let file = File::open(path)?;
        let mut rdr = csv::Reader::from_reader(file);

        let mut aln_vec: Vec<Aln> = rdr
            .deserialize::<Aln>()
            .filter_map(|res| res.ok()) // Skip invalid rows
            .filter(|aln| target_id.map_or(true, |tid| aln.target_id == tid))
            .collect();
        // sort by target_id and ref_start
        aln_vec.sort_by_key(|aln| (aln.target_id, aln.ref_start));

        Ok(Self::new(aln_vec.into_iter()))
    }
}

impl<I> Iterator for CsvPileupIterator<I>
where
    I: Iterator<Item = Aln>,
{
    type Item = (u32, u32, u32); // (target_id, position, depth)

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Load initial alignment if none loaded
            if self.current_target_id.is_none() {
                if let Some(aln) = self.input.peek().cloned() {
                    self.current_target_id = Some(aln.target_id);
                    self.pos = aln.ref_start;
                } else {
                    return None;
                }
            }

            let tid = self.current_target_id?;

            // Push new alignments into the active queue
            while let Some(aln) = self.input.peek() {
                if aln.target_id != tid || aln.ref_start > self.pos {
                    break;
                }
                self.active.push(self.input.next().unwrap());
            }

            // Remove alignments that have ended before current pos
            self.active.retain(|aln| self.pos <= aln.ref_end);

            // If no alignments and nothing left to load, finish
            if self.active.is_empty() {
                // Try to find the next target_id
                if let Some(aln) = self.input.peek() {
                    self.current_target_id = Some(aln.target_id);
                    self.pos = aln.ref_start;
                    continue;
                } else {
                    return None;
                }
            }

            // Yield current position's depth
            let result = Some((tid, self.pos, self.active.len() as u32));
            self.pos += 1;
            return result;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::create_dir_all;

    use crate::shared_test_functions::{compare_files, TEST_DATA, TEST_DIR};

    #[test]
    fn test_csv_pileup_iterator() {
        let alns_file = format!("{}/{}", TEST_DATA, "pileup/alns.csv");
        let expected_depths_file = format!("{}/{}", TEST_DATA, "pileup/expected_depths.csv");
        let iter = CsvPileupIterator::from_path(alns_file, None).expect("Failed to create iterator");

        create_dir_all(format!("{}/{}", TEST_DIR, "pileup"))
            .expect("Failed to create output directory");
        let output_file = format!("{}/{}", TEST_DIR, "pileup/depths.csv");
        let mut output =
            csv::Writer::from_path(&output_file).expect("Failed to create output file");
        for (target_id, pos, depth) in iter {
            output
                .serialize((target_id, pos, depth))
                .expect("Failed to write to output file");
        }

        output.flush().expect("Failed to flush output file");

        assert!(compare_files(&expected_depths_file, &output_file))
    }

    #[test]
    fn test_csv_pileup_iterator_target_id() {
        let alns_file = format!("{}/{}", TEST_DATA, "pileup/alns.csv");
        let iter = CsvPileupIterator::from_path(alns_file, Some(2))
            .expect("Failed to create iterator with target_id");

        let depths = iter.collect::<Vec<_>>();
        assert_eq!(depths.len(), 10);
        for (target_id, _pos, depth) in depths {
            assert_eq!(target_id, 2);
            assert_eq!(depth, 1);
        }
    }
}
