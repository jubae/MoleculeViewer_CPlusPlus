#include "MoleculeFactory.h"
#include "mainwindow.h"

#include <QDebug>

MoleculeFactory::MoleculeFactory(Ui::MainWindow& UI) : ui(&UI)
{
	qDebug() << "Molecule factory Initialised";
}

void MoleculeFactory::CreateMolecule(MoleculeData& NewMoleculeToCreate)
{
	qDebug() << NewMoleculeToCreate.MoleculeName;
}

