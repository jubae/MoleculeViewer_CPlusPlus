#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "moleculeglwidget.h"
#include "DataStructs.h"

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void onLoadMolecule();

private:
    Ui::MainWindow *ui;
    MoleculeData LoadedMolecule;
    MoleculeGLWidget* moleculeGLWidget; // Add this member
};
#endif // MAINWINDOW_H
