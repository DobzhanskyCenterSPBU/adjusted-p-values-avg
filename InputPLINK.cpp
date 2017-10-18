//
// Created by EL on 12.10.17.
//

#include "InputPLINK.h"
#include <plinkio/plinkio.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <fstream>

using namespace std;

InputPLINK::InputPLINK (char* fileName){
    this->filename = fileName;
}

// Move pointer to the appropriate row/SNP position (new_index) in .bed file
void bed_set_row(struct pio_file_t *plink_file, int new_index){
    size_t offset = bed_header_data_offset(&(&plink_file->bed_file)->header) + (new_index-1)*ceil(pio_num_samples(plink_file)/4);
    fseek((&plink_file->bed_file)->fp, offset, SEEK_SET);
    (&plink_file->bed_file)->cur_row = new_index;
}

vector<vector<unsigned short>> InputPLINK::createGenotypeMatrix(int lowInd, int upInd) {

    if (lowInd < 0) lowInd = 0;
    if (upInd < 0) upInd = 0;

    vector<vector<unsigned short>> result;
    struct pio_file_t plink_file; // Opens all three files at once (.bed, .bim, .fam)
    snp_t *snp_buffer;

    if (pio_open(&plink_file, this->filename) != PIO_OK) {
        cout << "Error here: Could not open " << this->filename << endl;
    }
    if(!pio_one_locus_per_row(&plink_file)) {
        cout << "This script requires that snps are rows and samples columns." << endl;
    }

    snp_buffer = (snp_t *) malloc( pio_row_size( &plink_file ));
    vector<unsigned short> snp_row;

    bed_set_row(&plink_file, lowInd); // Set the first SNP id we want to read

    int num_snps = 0; // Make sure we read the right amount of SNPs
    // Loop reads the interval [lowInd,upInd] (numeration starts at zero!) from genotype (.bed file)
    while(pio_next_row( &plink_file, snp_buffer ) == PIO_OK && num_snps < upInd - lowInd + 1){
        snp_row = {};
        for (int sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++){
            snp_row.push_back((unsigned short)snp_buffer[sample_id]);
        }
        result.push_back(snp_row);
        num_snps++;
    }

    free(snp_buffer);
    pio_close(&plink_file);

    return result;
}

vector<unsigned short> InputPLINK::createPhenotypeVector(){

    string fam_file_name(this->filename);
    fam_file_name += ".fam";

    vector<unsigned short> result;

    char delim = ' ';
    string fam_line_word;
    ifstream handle(fam_file_name);

    string fam_line;
    int ind;
    while (getline(handle, fam_line))
    {
        stringstream ss;
        ss.str(fam_line);
        ind = 0;
        while (getline(ss, fam_line_word, delim)) {
            ind++;
            if (ind == 6) {
                result.push_back((unsigned short) stoi(fam_line_word));
            }
        }
    }

    handle.close();

    return result;
}

/*
    struct pio_file_t plink_file;
    int sample_id;

    if( pio_open( &plink_file, this->filename ) != PIO_OK ) {
    }
    if( !pio_one_locus_per_row( &plink_file ) ) {
    }

    struct pio_sample_t *sample;
    for(sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++)
    {
        sample = pio_get_sample( &plink_file, sample_id );
        cout << sample->phenotype << " ";
    }
    pio_close( &plink_file );

 */