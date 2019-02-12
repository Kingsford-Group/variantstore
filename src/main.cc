/*
 * ============================================================================
 *
 *       Filename:  main.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */


#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <string>

#include "vcflib/Variant.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	vcflib::VariantCallFile variantFile;
	std::string filename = "./data/sample.vcf";
	variantFile.open(filename);
	vcflib::Variant var(variantFile);

	long int count = 0;
	while (variantFile.getNextVariant(var)) {
		count+= 1;
	}
	assert(count == 9);
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
