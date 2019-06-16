/*
 * ============================================================================
 *
 *       Filename:  bm_query.h
 *
 *         Author:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 										Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>

#include "variant_graph.h"
#include "dot_graph.h"
#include "query.h"

using namespace variantdb;
std::shared_ptr<spdlog::logger> console;

int main ( int argc, char *argv[] )
{
  console = spdlog::default_logger();

  std::string ref_file("/mnt/disk34/user/ppandey/vcfdata/ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.chromosome.22.fa");
	std::string vcf_file("/mnt/disk34/user/ppandey/vcfdata/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz");
	std::vector<std::string> vcfs = {vcf_file};
	PRINT("Creating variant graph");
	VariantGraph vg(ref_file, vcfs);

	PRINT("Creating Index");
	Index idx(&vg);
  /*
  // Read from serialized data
  PRINT("Reading chr22 ...");
  auto start = std::chrono::high_resolution_clock::now();
  VariantGraph vg("./chr22");
  PRINT("Reading chr22 ...");
  Index idx("./chrY");
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  PRINT("The time taken for loading is: " << duration.count()	<< "seconds");
  */
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // Take 10 random samples
  std::vector <string> samples;
  uint32_t max = 600;
  for (auto i=0; i<10; i++)
  {
    uint32_t id = rand() % max + 1;
    std::string sample_id = vg.get_sample_name(id);
    samples.push_back(sample_id);
  }

  // Generate 10 random positions
  std::vector <uint64_t> positions;
  uint64_t chr_len = vg.get_ref_length();
  PRINT("chr_len is: " << chr_len);
  uint64_t avg_query_len = chr_len / 10;
  for (auto i=0; i<10; i++)
  {
    uint64_t pos = rand() % chr_len + 1;
    positions.push_back(pos);
  }

  // Benchmark for query get_sample_var_in_ref
  PRINT("Benchmark for query get_sample_var_in_ref ...");
  start = std::chrono::high_resolution_clock::now();
  for (auto sample=samples.begin(); sample!=samples.end();++sample)
  {
    std::string sample_id = *sample;
    PRINT(sample_id);
    for (auto pos=positions.begin(); pos!=positions.end();++pos)
    {
      uint64_t pos_x = *pos;
      uint64_t pos_y = pos_x + avg_query_len;
      PRINT(pos_x << "," << pos_y);
      get_sample_var_in_ref(&vg, &idx, pos_x, pos_y, sample_id);
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  auto avg_time = duration.count() / 100;
  PRINT("The average time taken for get_sample_var_in_ref is: "
        << avg_time	<< "milliseconds");

  // Benchmark for query query_sample_from_ref
  PRINT("Benchmark for query_sample_from_ref ...");
  start = std::chrono::high_resolution_clock::now();
  for (auto sample=samples.begin(); sample!=samples.end();++sample)
  {
    std::string sample_id = *sample;
    PRINT(sample_id);
    for (auto pos=positions.begin(); pos!=positions.end();++pos)
    {
      uint64_t pos_x = *pos;
      uint64_t pos_y = pos_x + avg_query_len;
      PRINT(pos_x << "," << pos_y);
      query_sample_from_ref(&vg, &idx, pos_x, pos_y, sample_id);
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  avg_time = duration.count() / 100;
  PRINT("The average time taken for query_sample_from_ref is: "
        << avg_time	<< "milliseconds");

  // Benchmark for query query_sample_from_ref
  PRINT("Benchmark for query closest_var ...");
  start = std::chrono::high_resolution_clock::now();
  for (auto sample=samples.begin(); sample!=samples.end();++sample)
  {
    std::string sample_id = *sample;
    PRINT(sample_id);
    for (auto pos=positions.begin(); pos!=positions.end();++pos)
    {
      uint64_t pos_x = *pos;
      PRINT(pos_x);
      Variant var;
  		closest_var (&vg, &idx, pos_x, var);
    }
  }
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  avg_time = duration.count() / 100;
  PRINT("The average time taken for closest_var is: "
        << avg_time	<< "milliseconds");


  return 0;
}
