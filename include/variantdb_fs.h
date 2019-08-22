/*
 * ============================================================================
 *
 *       Filename:  variantdb_fs.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __VARIANTDB_FS_HPP__
#define __VARIANTDB_FS_HPP__

#include <vector>
#include <string>

namespace variantdb {
	namespace fs {
		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool FileExists(const char* path);
		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool DirExists(const char* path);
		void MakeDir(const char* path);
		// Taken from
		// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
		std::vector<std::string> GetFilesExt(const char *dir, const char *ext);
	}
}

#endif //__VARIANTDB_FS_HPP__
