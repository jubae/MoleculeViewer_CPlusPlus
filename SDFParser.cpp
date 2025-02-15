#include "SDFParser.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <map>

#include <QString>
#include <QDebug>


bool SDFParser::Read(const std::string& path, MoleculeData& NewMolecule)
{
	bool Status = false;

	std::ifstream stream(path);

	if (!stream.is_open())
	{
		qDebug() << "Unable to open file: " << path;
		return Status;
	}

	std::string line;
	int lineNumber = 0;

	// Skip Header Block
	for (int i = 0; i < 3; ++i)
	{
		std::getline(stream, line);
		lineNumber++;
	}

	// Table Reading
	int nAtoms = 0;
	int nBonds = 0;

	if (!std::getline(stream, line))
	{
		qDebug() << "Error reading counts line";
		stream.close();
		return Status;
	}

	lineNumber++;

	std::istringstream iss(line);

	if (!(iss >> nAtoms >> nBonds))
	{
		qDebug() << "Error parsing atom and bond counts from line " << lineNumber;
		stream.close();
		return Status;
	}

	qDebug() << "Number of Atoms: " << nAtoms;
	qDebug() << "Number of Bonds: " << nBonds;

	for (int i = 0; i < nAtoms; ++i)
	{
		std::getline(stream, line);
		lineNumber++;

		std::vector<std::string> AtomLine = SplitLine(line);

		if (!AtomLine.empty())
		{
			AtomData AtomDataToStore;

			AtomDataToStore.AtomNumber = i + 1;
			AtomDataToStore.AtomPosition = QVector3D(std::stof(AtomLine[0]), std::stof(AtomLine[1]), std::stof(AtomLine[2]));
			AtomDataToStore.AtomType = GetAtomType(AtomLine[3]);
			AtomDataToStore.AtomScale = 0.4f;

			NewMolecule.Atoms.push_back(AtomDataToStore);
		}
	}

	for (int i = 0; i < nBonds; ++i)
	{
		std::getline(stream, line);
		lineNumber++;
		std::vector<std::string> BondLine = SplitLine(line);

		if (!BondLine.empty())
		{
			BondData BondDataToStore;

			BondDataToStore.BondNumber = i + 1;
			BondDataToStore.Atom1Number = std::stof(BondLine[0]);
			BondDataToStore.Atom2Number = std::stof(BondLine[1]);
			BondDataToStore.BondType = GetBondType(std::stof(BondLine[2]));

			NewMolecule.Bonds.push_back(BondDataToStore);
		}
	}

	// Populate Molecule Name
	while (std::getline(stream, line))
	{
		std::string trimmedLine = line;
		trimmedLine.erase(0, trimmedLine.find_first_not_of(" \t"));
		trimmedLine.erase(trimmedLine.find_last_not_of(" \t") + 1);

		if (trimmedLine == "> <NAME>")
		{
			if (std::getline(stream, line))
			{
				std::string nameLine = line;
				nameLine.erase(0, nameLine.find_first_not_of(" \t"));
				nameLine.erase(nameLine.find_last_not_of(" \t") + 1);

				NewMolecule.MoleculeName = nameLine.c_str();
				break;
			}
		}
	}

	stream.close();
	Status = true;

	return Status;
}

std::vector<std::string> SDFParser::SplitLine(const std::string& line)
{
	std::vector<std::string> tokens;
	std::istringstream iss(line);
	std::string token;

	while (iss >> std::ws >> token)
	{
		tokens.push_back(token);
	}

	return tokens;
}

EAtomType SDFParser::GetAtomType(const std::string& AtomSymbol)
{
	// remove case sensitivity
	std::string lowerSymbol = AtomSymbol;
	std::transform(lowerSymbol.begin(), lowerSymbol.end(), lowerSymbol.begin(), ::tolower); 

	std::map<std::string, EAtomType> atomTypeMap =
	{
		{"c", EAtomType::Carbon},
		{"h", EAtomType::Hydrogen},
		{"n", EAtomType::Nitrogen},
		{"o", EAtomType::Oxygen},
		{"f", EAtomType::Fluorine},
		{"si", EAtomType::Silicon},
		{"p", EAtomType::Phosphorus},
		{"s", EAtomType::Sulphur},
		{"cl", EAtomType::Chlorine},
		{"cu", EAtomType::Copper},
		{"zn", EAtomType::Zinc},
		{"br", EAtomType::Bromine},
		{"ag", EAtomType::Silver},
		{"i", EAtomType::Iodine},
		{"au", EAtomType::Gold},
		{"hg", EAtomType::Mercury},
	};

	auto it = atomTypeMap.find(lowerSymbol);

	if (it != atomTypeMap.end())
	{
		return it->second;
	}
	else
	{
		return EAtomType::Unknown;
	}
}

EBondType SDFParser::GetBondType(const int& BondSymbol)
{
	std::map<int, EBondType> BondTypeMap =
	{
		{1, EBondType::Single},
		{2, EBondType::Double},
		{3, EBondType::Triple}
	};

	auto it = BondTypeMap.find(BondSymbol);

	if (it != BondTypeMap.end())
	{
		return it->second;
	}
	else
	{
		return EBondType::Unknown;
	}
}

std::string SDFParser::getMoleculeName(const std::string& path)
{
	std::ifstream stream(path);
	std::string moleculeName = "";

	if (stream.is_open())
	{
		std::string line;

		while (std::getline(stream, line))
		{
			std::string trimmedLine = line;
			trimmedLine.erase(0, trimmedLine.find_first_not_of(" \t"));
			trimmedLine.erase(trimmedLine.find_last_not_of(" \t") + 1);

			if (trimmedLine == "> <NAME>")
			{
				if (std::getline(stream, line))
				{
					std::string nameLine = line;
					nameLine.erase(0, nameLine.find_first_not_of(" \t"));
					nameLine.erase(nameLine.find_last_not_of(" \t") + 1);

					moleculeName = nameLine;
					break;
				}
			}
		}
		stream.close();
	}
	else
	{
		std::cerr << "Unable to open file: " << path << std::endl;
	}

	return moleculeName;
}