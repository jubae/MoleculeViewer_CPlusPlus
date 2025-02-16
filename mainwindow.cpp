#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "SDFParser.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QApplication>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->actionLoadMolecule, &QAction::triggered, this, &MainWindow::onLoadMolecule);
    connect(ui->actionExit, &QAction::triggered, this, &MainWindow::onExit);
    connect(ui->actionAbout, &QAction::triggered, this, &MainWindow::about);

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
            ui->MoleculeWidget->loadMolecule(LoadedMolecule);
        }
        else
        {
            QMessageBox::information(nullptr, "Information", "Unable to parse the selected file.");
        }

    }
}

void MainWindow::onExit()
{
    QApplication::quit();
}

void MainWindow::about()
{
    QMessageBox msgBox;
    msgBox.setText("<div align=\"center\"><b>Molecule Viewer V1.0.0</b><br>Created as a demonstration by Mo <b>Davis-Afsar </b>2025<br>Created Using C++ and Qt</div>");
    msgBox.setWindowTitle("Molecule Viewer");
    msgBox.setWindowFlags(msgBox.windowFlags() & ~Qt::WindowContextHelpButtonHint);
    msgBox.setIcon(QMessageBox::NoIcon);
    msgBox.exec();
}

