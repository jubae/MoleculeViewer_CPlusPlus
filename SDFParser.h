#ifndef SDFPARSER_H
#define SDFPARSER_H

#include "DataStructs.h"

#include <string>
#include <vector>

class SDFParser
{
public:
	bool Read(const std::string& path, MoleculeData& NewMolecule);

private:
	std::vector<std::string> SplitLine(const std::string& line);
	EAtomType GetAtomType(const std::string& AtomSymbol);
	EBondType GetBondType(const int& BondSymbol);
	std::string getMoleculeName(const std::string& path);
};

#endif // SDFPARSER_H

