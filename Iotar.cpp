// Copyright 2020      Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Iotar.h"

#include <sstream>
#include <archive.h>
#include <archive_entry.h>

using namespace std;

// Mixture of C (for libarchive code) and C++ (for the interface)

map<string, istream*> ifstream_from_tar(string tarfile) {
	map<string, istream*> mymap;
	
	struct archive *a;
	struct archive_entry *entry;
	int r;

	a = archive_read_new();
	//archive_read_support_compression_all(a);
	archive_read_support_filter_all(a);
	archive_read_support_format_all(a);
	
	r = archive_read_open_filename(a, tarfile.c_str(), 1024*1024); 
	if (r != ARCHIVE_OK) {
		cerr << "Problem opening the archive " << tarfile << endl;
		exit(EXIT_FAILURE);
	}
	while (archive_read_next_header(a, &entry) == ARCHIVE_OK) {
		long int entry_size = archive_entry_size(entry);
		char *fileContents = (char*)malloc(entry_size+1);
		fileContents[entry_size] = '\0'; // if not, istringstream reads beyond!!!
		archive_read_data(a, fileContents, entry_size);
		
		string filename (archive_entry_pathname(entry));
		
		mymap.insert(pair<string, istream*>(filename,  new istringstream(fileContents)));
		//~ cerr << filename << "(" << entry_size << ") :" << endl;
		//~ cerr << "*************************" << endl;
		//~ cerr << fileContents << endl;
		//~ cerr << "**************************" << endl;
		//~ archive_entry_clear(entry);
		//~ archive_read_data_skip(a);
		free(fileContents); // because stringstream copies its own internal buffer
	}
	r = archive_read_free(a);
	if (r != ARCHIVE_OK) {
		cerr << "Problem closing the archive " << tarfile << endl;
		exit(EXIT_FAILURE);
	}
	return mymap;
}

// This is a helper function
string my_basename(string filePath, char seperator)
{
    size_t sepPos = filePath.rfind(seperator);
    if(sepPos != string::npos)
    {
        return filePath.substr(sepPos + 1, filePath.size());
    }
    return filePath;
}
