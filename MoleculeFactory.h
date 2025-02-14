#ifndef MOLECULEFACTORY_H
#define MOLECULEFACTORY_H

#include "MainWindow.h"
#include "DataStructs.h"

class MoleculeFactory
{
public:
	MoleculeFactory(Ui::MainWindow& UI);

	void CreateMolecule(MoleculeData& NewMoleculeToCreate);

private:   
	Ui::MainWindow* ui;
};

#endif //MOLECULEFACTORY_H

