#include "moleculeglwidget.h"

#include <algorithm>
#include <limits>

#include <QWheelEvent>
#include <QMouseEvent>
#include <QDebug>
#include <QtMath>
#include <QPainter>


MoleculeGLWidget::MoleculeGLWidget(QWidget* parent) : QOpenGLWidget(parent)
{
    setFocusPolicy(Qt::StrongFocus);
}

MoleculeGLWidget::~MoleculeGLWidget()
{
    makeCurrent();
    delete program;
    if (atomVao && atomVao->isCreated())
    {
        atomVao->destroy();
        delete atomVao;
    }
    if (atomVbo) glDeleteBuffers(1, &atomVbo);
    if (atomEbo) glDeleteBuffers(1, &atomEbo);

    if (bondVao && bondVao->isCreated()) 
    {
        bondVao->destroy();
        delete bondVao;
    }
    if (bondVbo) glDeleteBuffers(1, &bondVbo);
    if (bondEbo) glDeleteBuffers(1, &bondEbo);


    doneCurrent();
}

void MoleculeGLWidget::loadMolecule(const MoleculeData& MoleculeToGenerate)
{
    moleculeData = MoleculeToGenerate;
    atomMeshes.clear();
    for (const auto& atom : moleculeData.Atoms) 
    {
        atomMeshes.append(createAtom(atom.AtomScale, 32, 32));
    }

    centerAtoms();

    cylinderMesh = createCylinder(0.025f, 1.0f, 16);

    update();
}


AtomMeshData MoleculeGLWidget::createAtom(float radius, int slices, int stacks)
{
    AtomMeshData meshData;
    for (int i = 0; i <= stacks; ++i) {
        float phi = M_PI * i / static_cast<float>(stacks);
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * M_PI * j / static_cast<float>(slices);
            float x = radius * std::sin(phi) * std::cos(theta);
            float y = radius * std::cos(phi);
            float z = radius * std::sin(phi) * std::sin(theta);
            meshData.vertices.append(x);
            meshData.vertices.append(y);
            meshData.vertices.append(z);
            QVector3D normal(x, y, z);
            normal.normalize();
            meshData.normals.append(normal.x());
            meshData.normals.append(normal.y());
            meshData.normals.append(normal.z());
        }
    }
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            GLuint index1 = i * (slices + 1) + j;
            GLuint index2 = index1 + 1;
            GLuint index3 = index1 + slices + 1;
            GLuint index4 = index3 + 1;
            meshData.indices.append(index1);
            meshData.indices.append(index2);
            meshData.indices.append(index3);
            meshData.indices.append(index2);
            meshData.indices.append(index4);
            meshData.indices.append(index3);
        }
    }
    return meshData;
}

AtomMeshData MoleculeGLWidget::createCylinder(float radius, float height, int slices)
{
    AtomMeshData meshData;

    for (int i = 0; i <= slices; ++i) 
    {
        float theta = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(slices);
        float x = radius * std::cos(theta);
        float z = radius * std::sin(theta);

        // --- Bottom vertex ---
        meshData.vertices.append(x);
        meshData.vertices.append(-height / 2.0f);
        meshData.vertices.append(z);
        QVector3D normal(x, 0, z);
        normal.normalize();
        meshData.normals.append(normal.x());
        meshData.normals.append(normal.y());
        meshData.normals.append(normal.z());

        // --- Top vertex ---
        meshData.vertices.append(x);
        meshData.vertices.append(height / 2.0f);
        meshData.vertices.append(z);
        normal = QVector3D(x, 0, z);
        normal.normalize();
        meshData.normals.append(normal.x());
        meshData.normals.append(normal.y());
        meshData.normals.append(normal.z());
    }


    for (int i = 0; i < slices; ++i) 
    {
        GLuint index1 = 2 * i;
        GLuint index2 = 2 * i + 1;
        GLuint index3 = 2 * (i + 1); 
        GLuint index4 = 2 * (i + 1) + 1;


        if (i == slices - 1)
        {
            index3 = 0;
            index4 = 1;
        }

        meshData.indices.append(index1);
        meshData.indices.append(index3);
        meshData.indices.append(index2);

        meshData.indices.append(index2);
        meshData.indices.append(index3);
        meshData.indices.append(index4);
    }

    return meshData;
}


void MoleculeGLWidget::initializeGL()
{
    initializeOpenGLFunctions();
    program = new QOpenGLShaderProgram(this);

    atomColourMap[EAtomColour::White] = QVector3D(1.0f, 1.0f, 1.0f);
    atomColourMap[EAtomColour::Red] = QVector3D(1.0f, 0.0f, 0.0f);
    atomColourMap[EAtomColour::Green] = QVector3D(0.0f, 1.0f, 0.0f);
    atomColourMap[EAtomColour::Blue] = QVector3D(0.0f, 0.0f, 1.0f);
    atomColourMap[EAtomColour::Grey] = QVector3D(0.5f, 0.5f, 0.5f);
    atomColourMap[EAtomColour::Yellow] = QVector3D(1.0f, 1.0f, 0.0f);
    atomColourMap[EAtomColour::Orange] = QVector3D(1.0f, 0.5f, 0.0f);
    atomColourMap[EAtomColour::DarkRed] = QVector3D(0.5f, 0.0f, 0.0f);

    // --- Vertex Shader ---
    program->addShaderFromSourceCode(QOpenGLShader::Vertex,
        "#version 330 core\n"
        "layout (location = 0) in vec3 position;\n"
        "layout (location = 1) in vec3 normal;\n"
        "uniform mat4 matrix;\n"
        "uniform mat4 projection;\n"
        "out vec3 fragNormal;\n"
        "out vec3 fragPosWorld;\n"
        "out vec3 fragPosLocal;\n"
        "void main() {\n"
        "   gl_Position = projection * matrix * vec4(position, 1.0);\n"
        "   fragNormal = mat3(transpose(inverse(matrix))) * normal;\n"
        "   fragPosWorld = vec3(matrix * vec4(position, 1.0));\n"
        "   fragPosLocal = position;\n"
        "}\n");

    // --- Fragment Shader ---
    program->addShaderFromSourceCode(QOpenGLShader::Fragment,
        "#version 330 core\n"
        "in vec3 fragNormal;\n"
        "in vec3 fragPosWorld;\n"
        "in vec3 fragPosLocal;\n"
        "uniform vec3 atomColour;\n"
        "uniform vec3 atomColour2;\n"
        "uniform vec3 lightDir;\n"
        "uniform vec3 lightColour;\n"
        "uniform vec3 viewPos;\n"
        "out vec4 fragColour;\n"
        "void main() {\n"
        "   float ambientStrength = 0.3;\n"
        "   vec3 ambient = ambientStrength * lightColour;\n"
        "   vec3 norm = normalize(fragNormal);\n"
        "   float diff = max(dot(norm, lightDir), 0.0);\n"
        "   vec3 diffuse = diff * lightColour;\n"
        "   float specularStrength = 0.5;\n"
        "   vec3 viewDir = normalize(viewPos - fragPosWorld);\n"
        "   vec3 reflectDir = reflect(-lightDir, norm);\n"
        "   float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);\n"
        "   vec3 specular = specularStrength * spec * lightColour;\n"
        "   vec3 result = (ambient + diffuse + specular) * (fragPosLocal.y < 0.0 ? atomColour : atomColour2);\n"
        "   fragColour = vec4(result, 1.0);\n"
        "}\n");

    if (!program->link())
    {
        qDebug() << "Shader link error:" << program->log();
        return;
    }

    atomVao = new QOpenGLVertexArrayObject(this);
    atomVao->create();
    atomVao->bind();
    glGenBuffers(1, &atomVbo);
    glBindBuffer(GL_ARRAY_BUFFER, atomVbo);

    program->enableAttributeArray(0); // Position
    program->setAttributeBuffer(0, GL_FLOAT, 0, 3, 6 * sizeof(GLfloat));
    program->enableAttributeArray(1); // Normal
    program->setAttributeBuffer(1, GL_FLOAT, 3 * sizeof(GLfloat), 3, 6 * sizeof(GLfloat));

    glGenBuffers(1, &atomEbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, atomEbo);
    atomVao->release();

    // --- Bond VAO Setup ---
    bondVao = new QOpenGLVertexArrayObject(this);
    bondVao->create();
    bondVao->bind();
    glGenBuffers(1, &bondVbo);
    glBindBuffer(GL_ARRAY_BUFFER, bondVbo);

    program->enableAttributeArray(0); // Position
    program->setAttributeBuffer(0, GL_FLOAT, 0, 3, 6 * sizeof(GLfloat));
    program->enableAttributeArray(1); // Normal
    program->setAttributeBuffer(1, GL_FLOAT, 3 * sizeof(GLfloat), 3, 6 * sizeof(GLfloat));

    glGenBuffers(1, &bondEbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bondEbo);
    bondVao->release();

    moleculeMatrix.setToIdentity();
    modelMatrix.setToIdentity();
    centeringTranslation = QVector3D(0, 0, 0);

    QMatrix4x4 projection;
    projection.perspective(45.0f, 1.0f, 0.1f, 100.0f);
    program->bind();
    program->setUniformValue("projection", projection);
    QVector3D lightDir = QVector3D(1.0f, 1.0f, 1.0f).normalized();
    program->setUniformValue("lightDir", lightDir);

    // --- Add light Colour ---
    QVector3D lightColour(1.0f, 1.0f, 1.0f); // White light
    program->setUniformValue("lightColour", lightColour);

    program->release();

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
}

void MoleculeGLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    program->bind();
    QVector3D lightDir = QVector3D(1.0f, 1.0f, 1.0f).normalized();
    program->setUniformValue("lightDir", lightDir);

    // --- Camera position to shader ---
    QVector3D viewPos = modelMatrix.inverted().map(QVector3D(0.0f, 0.0f, 0.0f));
    program->setUniformValue("viewPos", viewPos);

    // --- Render Atoms ---
    if (!moleculeData.Atoms.isEmpty() && atomMeshes.size() == moleculeData.Atoms.size())
    {
        for (int i = 0; i < moleculeData.Atoms.size(); ++i)
        {
            const auto& atom = moleculeData.Atoms[i];
            const auto& atomMesh = atomMeshes[i];
            EAtomColour atomColour = getAtomColour(atom.AtomType);
            QVector3D Colour = atomColourMap[atomColour];
            program->setUniformValue("atomColour", Colour);
            program->setUniformValue("atomColour2", Colour);

            // --- Atom scale ---
            QMatrix4x4 finalMatrix = moleculeMatrix * modelMatrix;
            finalMatrix.translate(atom.AtomPosition);
            finalMatrix.scale(atom.AtomScale);
            program->setUniformValue("matrix", finalMatrix);

            atomVao->bind();
            glBindBuffer(GL_ARRAY_BUFFER, atomVbo);

            QVector<GLfloat> combinedData;
            for (int j = 0; j < atomMesh.vertices.size() / 3; ++j)
            {
                // Vertex Position
                combinedData.append(atomMesh.vertices[j * 3]);
                combinedData.append(atomMesh.vertices[j * 3 + 1]);
                combinedData.append(atomMesh.vertices[j * 3 + 2]);

                // Normal Vector
                combinedData.append(atomMesh.normals[j * 3]);
                combinedData.append(atomMesh.normals[j * 3 + 1]);
                combinedData.append(atomMesh.normals[j * 3 + 2]);
            }

            glBufferData(GL_ARRAY_BUFFER, combinedData.size() * sizeof(GLfloat), combinedData.data(), GL_STATIC_DRAW);

            program->enableAttributeArray(0);
            program->setAttributeBuffer(0, GL_FLOAT, 0, 3, 6 * sizeof(GLfloat));
            program->enableAttributeArray(1);
            program->setAttributeBuffer(1, GL_FLOAT, 3 * sizeof(GLfloat), 3, 6 * sizeof(GLfloat));

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, atomEbo);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, atomMesh.indices.size() * sizeof(GLuint), atomMesh.indices.data(), GL_STATIC_DRAW);
            glDrawElements(GL_TRIANGLES, atomMesh.indices.size(), GL_UNSIGNED_INT, 0);
            atomVao->release();
        }
    }

    renderBonds();

    program->release();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width(), height(), 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (moleculeData.Atoms.size()) // --- If a Molecule is loaded ---
    {
        // Molecule Name
        QPainter painter;
        painter.begin(this);
        painter.setPen(Qt::white);
        painter.setFont(QFont("Arial", 24));

        QRect textRect(20, 20, width() - 40, 35);
        painter.drawText(textRect, Qt::AlignLeft | Qt::AlignTop, moleculeData.MoleculeName);

        // Instructions
        painter.setFont(QFont("Arial", 12));
        painter.setPen(Qt::lightGray);

        QRect LMBtextRect(20, 100, width() - 40, 30);
        painter.drawText(LMBtextRect, Qt::AlignLeft | Qt::AlignTop, "Left Click : Rotate Molecule");


        QRect MMBtextRect(20, 120, width() - 40, 30);
        painter.drawText(MMBtextRect, Qt::AlignLeft | Qt::AlignTop, "Middle click : Pan Molecule");


        QRect RMBtextRect(20, 140, width() - 40, 30);
        painter.drawText(RMBtextRect, Qt::AlignLeft | Qt::AlignTop, "Right click : Reset Molecule");

        painter.setFont(QFont("Arial", 12));

        QRect WheeltextRect(20, 160, width() - 40, 30);
        painter.drawText(WheeltextRect, Qt::AlignLeft | Qt::AlignTop, "Mouse Wheel : Zoom Molecule");

        // Atoms Number
        painter.setPen(Qt::green);
        painter.setFont(QFont("Courier New", 10, QFont::Bold));
        QRect anotherTextRect(10, 50, width() - 20, 20);
        painter.drawText(anotherTextRect, Qt::AlignRight | Qt::AlignTop, QString("Atoms: %1").arg(moleculeData.Atoms.size()));

        painter.end();
    }
    else
    {
        QPainter painter;
        painter.begin(this);
        painter.setPen(Qt::white);
        painter.setFont(QFont("Courier New", 24, QFont::Bold));

        // Main Load Molecule Message
        QRect textRect(0, 0, width(), height());
        painter.drawText(textRect, Qt::AlignCenter, "Use the File menu to \n load a molecule");

        painter.end();
    }

    glPopMatrix();
    glPopAttrib();

    GLenum error = glGetError();
    if (error != GL_NO_ERROR)
    {
        qDebug() << "OpenGL error:" << error;
    }
}

void MoleculeGLWidget::renderBonds()
{
    bondVao->bind();
    glBindBuffer(GL_ARRAY_BUFFER, bondVbo);

    QVector<GLfloat> combinedData;
    for (int j = 0; j < cylinderMesh.vertices.size() / 3; ++j)
    {
        combinedData.append(cylinderMesh.vertices[j * 3]);
        combinedData.append(cylinderMesh.vertices[j * 3 + 1]);
        combinedData.append(cylinderMesh.vertices[j * 3 + 2]);
        combinedData.append(cylinderMesh.normals[j * 3]);
        combinedData.append(cylinderMesh.normals[j * 3 + 1]);
        combinedData.append(cylinderMesh.normals[j * 3 + 2]);
    }

    glBufferData(GL_ARRAY_BUFFER, combinedData.size() * sizeof(GLfloat), combinedData.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bondEbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, cylinderMesh.indices.size() * sizeof(GLuint), cylinderMesh.indices.data(), GL_STATIC_DRAW);

    for (const auto& bond : moleculeData.Bonds)
    {
        // --- Get atom indices from atom numbers ---
        int atomIndex1 = getAtomIndexFromNumber(bond.Atom1Number);
        int atomIndex2 = getAtomIndexFromNumber(bond.Atom2Number);

        Q_ASSERT(atomIndex1 != -1);
        Q_ASSERT(atomIndex2 != -1);

        QVector3D pos1 = moleculeData.Atoms[atomIndex1].AtomPosition;
        QVector3D pos2 = moleculeData.Atoms[atomIndex2].AtomPosition;

        QVector3D Colour1 = atomColourMap[getAtomColour(moleculeData.Atoms[atomIndex1].AtomType)];
        QVector3D Colour2 = atomColourMap[getAtomColour(moleculeData.Atoms[atomIndex2].AtomType)];

        QVector3D bondVector = pos2 - pos1;
        float bondLength = bondVector.length();
        QVector3D bondDirection = bondVector.normalized();

        // --- Calculate Bond midpoint ---
        QVector3D midPoint = (pos1 + pos2) / 2.0f;

        QMatrix4x4 bondMatrix = moleculeMatrix * modelMatrix;

        bondMatrix.translate(midPoint);

        QQuaternion rotation = QQuaternion::rotationTo(QVector3D(0, 1, 0), bondDirection);
        bondMatrix.rotate(rotation);

        bondMatrix.scale(1.0f, bondLength, 1.0f); // Scale AFTER translation and rotation

        program->setUniformValue("matrix", bondMatrix);
        program->setUniformValue("atomColour", Colour1);
        program->setUniformValue("atomColour2", Colour2);

        int numBonds = 1;
        if (bond.BondType == EBondType::Double) numBonds = 2;
        if (bond.BondType == EBondType::Triple) numBonds = 3;

        float bondSpacing = 0.1f;

        for (int i = 0; i < numBonds; ++i)
        {
            QMatrix4x4 instanceMatrix = bondMatrix;

            if (numBonds > 1)
            {
                QVector3D perpendicular = QVector3D::crossProduct(bondDirection, QVector3D(0, 1, 0));

                if (perpendicular.lengthSquared() < 0.0001f)
                {
                    perpendicular = QVector3D::crossProduct(bondDirection, QVector3D(1, 0, 0));
                }
                perpendicular.normalize();

                float offset = (i - (numBonds - 1) / 2.0f) * bondSpacing;
                instanceMatrix.translate(perpendicular * offset);
            }

            program->setUniformValue("matrix", instanceMatrix);

            glDrawElements(GL_TRIANGLES, cylinderMesh.indices.size(), GL_UNSIGNED_INT, 0);
        }
    }

    bondVao->release();
}


void MoleculeGLWidget::resizeGL(int w, int h)
{
    float aspectRatio = static_cast<float>(w) / static_cast<float>(h ? h : 1);

    int side = qMin(w, h);

    int viewportX = (w - side) / 2;
    int viewportY = (h - side) / 2;

    glViewport(viewportX, viewportY, side, side);

    QMatrix4x4 projection;
    projection.perspective(45.0f, aspectRatio, 0.1f, 100.0f);

    program->bind();
    program->setUniformValue("projection", projection);
    program->release();
}

void MoleculeGLWidget::centerAtoms()
{
    if (moleculeData.Atoms.isEmpty()) return;

    QVector3D minPos(std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max());
    QVector3D maxPos(std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest());

    for (const auto& atom : moleculeData.Atoms)
    {
        minPos.setX(qMin(minPos.x(), atom.AtomPosition.x()));
        minPos.setY(qMin(minPos.y(), atom.AtomPosition.y()));
        minPos.setZ(qMin(minPos.z(), atom.AtomPosition.z()));
        maxPos.setX(qMax(maxPos.x(), atom.AtomPosition.x()));
        maxPos.setY(qMax(maxPos.y(), atom.AtomPosition.y()));
        maxPos.setZ(qMax(maxPos.z(), atom.AtomPosition.z()));
    }

    QVector3D center = (minPos + maxPos) / 2.0f;
    float maxExtent = qMax(maxPos.x() - minPos.x(), qMax(maxPos.y() - minPos.y(), maxPos.z() - minPos.z()));


    float scaleFactor = 1.0f / maxExtent;

    moleculeMatrix.setToIdentity();
    modelMatrix.setToIdentity();

    centeringTranslation = center;

    modelMatrix.translate(0.0f, 0.0f, -1.5f);  // Move camera BACK *first*
    modelMatrix.scale(scaleFactor);            // *Then* scale
    modelMatrix.translate(-center);            // *Then* center to get view to be right at start

    update();
}

void MoleculeGLWidget::wheelEvent(QWheelEvent* event)
{
    float zoomFactor = 1.0f + event->angleDelta().y() / 1000.0f;

    modelMatrix.translate(centeringTranslation);
    modelMatrix.scale(zoomFactor);
    modelMatrix.translate(-centeringTranslation);

    update();

    event->accept();
}

void MoleculeGLWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton || event->button() == Qt::MiddleButton)
    {
        lastMousePos = event->pos();
        event->accept();
    }
    else if (event->button() == Qt::RightButton)
    {
        centerAtoms();
        event->accept();
    }
}

void MoleculeGLWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (event->buttons() & Qt::MiddleButton)
    {
        QPointF delta = event->pos() - lastMousePos;
        lastMousePos = event->pos();

        float panScale = 0.005f;
        float panX = delta.x() * panScale;
        float panY = -delta.y() * panScale;

        modelMatrix.translate(panX, panY, 0.0f);

        update();
        event->accept();
    }

    if (event->buttons() & Qt::LeftButton)
    {
        QPointF delta = event->pos() - lastMousePos;
        lastMousePos = event->pos();

        float rotationScale = 0.5f;
        float rotationX = delta.y() * rotationScale;
        float rotationY = delta.x() * rotationScale;

        QQuaternion xRotation = QQuaternion::fromAxisAndAngle(QVector3D(1, 0, 0), rotationX);
        QQuaternion yRotation = QQuaternion::fromAxisAndAngle(QVector3D(0, 1, 0), rotationY);
        QQuaternion totalRotation = xRotation * yRotation;

        modelMatrix.translate(centeringTranslation);
        modelMatrix.rotate(totalRotation);
        modelMatrix.translate(-centeringTranslation);

        update();
        event->accept();
    }
}

EAtomColour MoleculeGLWidget::getAtomColour(EAtomType atomType)
{
    switch (atomType)
    {
        case EAtomType::Hydrogen: return EAtomColour::White;
        case EAtomType::Boron: return EAtomColour::DarkRed;
        case EAtomType::Carbon: return EAtomColour::Grey;
        case EAtomType::Nitrogen: return EAtomColour::Blue;
        case EAtomType::Oxygen: return EAtomColour::Red;
        case EAtomType::Fluorine: return EAtomColour::DarkRed;
        case EAtomType::Silicon: return EAtomColour::DarkRed;
        case EAtomType::Phosphorus: return EAtomColour::Orange;
        case EAtomType::Sulphur: return EAtomColour::Yellow;
        case EAtomType::Chlorine: return EAtomColour::Green;
        case EAtomType::Copper: return EAtomColour::DarkRed;
        case EAtomType::Zinc: return EAtomColour::DarkRed;
        case EAtomType::Bromine: return EAtomColour::DarkRed;
        case EAtomType::Silver: return EAtomColour::DarkRed;
        case EAtomType::Iodine: return EAtomColour::DarkRed;
        case EAtomType::Gold: return EAtomColour::DarkRed;
        case EAtomType::Mercury: return EAtomColour::DarkRed;
        case EAtomType::Metal: return EAtomColour::DarkRed;
        default: return EAtomColour::DarkRed; // if AtomType is unknown this is default colour
    }
}

QString MoleculeGLWidget::atomColourToString(EAtomColour colour)
{
    switch (colour)
    {
        case EAtomColour::White: return "White";
        case EAtomColour::Red: return "Red";
        case EAtomColour::Green: return "Green";
        case EAtomColour::Blue: return "Blue";
        case EAtomColour::Grey: return "Grey";
        case EAtomColour::Yellow: return "Yellow";
        case EAtomColour::Orange: return "Orange";
        case EAtomColour::DarkRed: return "DarkRed";
        default: return "Unknown";
    }
}

int MoleculeGLWidget::getAtomIndexFromNumber(int atomNumber) const
{
    for (int i = 0; i < moleculeData.Atoms.size(); ++i)
    {
        if (moleculeData.Atoms[i].AtomNumber == atomNumber)
        {
            return i;
        }
    }

    return -1;
}
