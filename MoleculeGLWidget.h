// MoleculeGLWidget.h
#ifndef MOLECULEGLWIDGET_H
#define MOLECULEGLWIDGET_H

#include "datastructs.h"

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QMatrix4x4>
#include <QVector3D>
#include <QMap>
#include <QPointF>
#include <QQuaternion>
#include <QPainter>

class MoleculeGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    MoleculeGLWidget(QWidget* parent = nullptr);
    ~MoleculeGLWidget();

    void loadMolecule(const MoleculeData& MoleculeToGenerate);

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;

    void wheelEvent(QWheelEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;

private:
    MoleculeData moleculeData;
    QOpenGLShaderProgram* program;

    AtomMeshData createAtom(float radius, int slices, int stacks);

    QOpenGLVertexArrayObject* atomVao;
    GLuint atomVbo;
    GLuint atomEbo;

    QOpenGLVertexArrayObject* bondVao;
    GLuint bondVbo;
    GLuint bondEbo;
    AtomMeshData cylinderMesh;

    QMatrix4x4 moleculeMatrix;
    QMatrix4x4 modelMatrix;
    QVector<AtomMeshData> atomMeshes;
    QPointF lastMousePos;
    QVector3D centeringTranslation;
    void centerAtoms();
    QMap<EAtomColour, QVector3D> atomColorMap;
    EAtomColour getAtomColor(EAtomType atomType);
    QString atomColourToString(EAtomColour colour);
    AtomMeshData createCylinder(float radius, float height, int slices);
    void renderBonds();
    int getAtomIndexFromNumber(int atomNumber) const; // Add this line
    int SelectedAtom = 0;
};

#endif // MOLECULEGLWIDGET_H