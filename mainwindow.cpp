#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "SDFParser.h"

#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    moleculeGLWidget = new MoleculeGLWidget(this); // Create the OpenGL widget
    setCentralWidget(moleculeGLWidget); // Set it as the central widget

    connect(ui->actionLoadMolecule, &QAction::triggered, this, &MainWindow::onLoadMolecule);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::onLoadMolecule()
{
    QString MoleculeFilePath = QFileDialog::getOpenFileName(this, "Load Molecule", "", "SDF Files (*.sdf)");
    if (!MoleculeFilePath.isEmpty())
    {
        LoadedMolecule.Atoms.clear();
        LoadedMolecule.Bonds.clear();
        LoadedMolecule.MoleculeName = "";

        SDFParser Parser;
        bool ReadSuccessful = Parser.Read(MoleculeFilePath.toStdString(), LoadedMolecule);

        if (ReadSuccessful)
        {
            moleculeGLWidget->loadMolecule(LoadedMolecule);
        }

    }
}
