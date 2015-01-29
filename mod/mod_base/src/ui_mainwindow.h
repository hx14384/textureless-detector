/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Mon Jan 24 09:48:45 2011
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSlider>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTabWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "GLWidget.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    GLWidget *mpGLWidget;
    QTabWidget *mpTabWidget;
    QWidget *Control;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout_5;
    QLabel *label;
    QPushButton *mpNextButton;
    QLabel *label2;
    QPushButton *mpNextObjButton;
    QPushButton *mpPrevObjButton;
    QCheckBox *mpIsLearnCheckBox;
    QSpacerItem *verticalSpacer;
    QWidget *Flags;
    QGroupBox *groupBox_4;
    QWidget *layoutWidget_3;
    QVBoxLayout *verticalLayout;
    QCheckBox *mpShowMessagesCheckBox;
    QCheckBox *mpShowEdgesCheckBox;
    QCheckBox *mpShowChainsCheckBox;
    QCheckBox *mpShowAllChainsCheckBox;
    QCheckBox *mpShowLinesCheckBox;
    QGroupBox *groupBox_5;
    QWidget *layoutWidget_4;
    QVBoxLayout *verticalLayout_3;
    QSlider *fps_slider; 

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
//        MainWindow->resize(1440, 960);
        MainWindow->resize(800, 480);
        QSizePolicy sizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
//        MainWindow->setMinimumSize(QSize(1440, 960));
        MainWindow->setMinimumSize(QSize(800, 480));
        MainWindow->setDockNestingEnabled(false);
        MainWindow->setDockOptions(QMainWindow::AnimatedDocks|QMainWindow::VerticalTabs);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(centralWidget->sizePolicy().hasHeightForWidth());
        centralWidget->setSizePolicy(sizePolicy1);
        mpGLWidget = new GLWidget(centralWidget);
        mpGLWidget->setObjectName(QString::fromUtf8("mpGLWidget"));
      //  mpGLWidget->setGeometry(QRect(0, 0, 1280, 960));
        mpGLWidget->setGeometry(QRect(0, 0, 640, 480));
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(mpGLWidget->sizePolicy().hasHeightForWidth());
        mpGLWidget->setSizePolicy(sizePolicy2);
        mpGLWidget->setMinimumSize(QSize(320, 240));
        mpTabWidget = new QTabWidget(centralWidget);
        mpTabWidget->setObjectName(QString::fromUtf8("mpTabWidget"));
        mpTabWidget->setGeometry(QRect(640, 0, 160, 480));
    //    mpTabWidget->setGeometry(QRect(1280, 0, 160, 960));
        Control = new QWidget();
        Control->setObjectName(QString::fromUtf8("Control"));
        layoutWidget = new QWidget(Control);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 0, 141, 441));
        verticalLayout_5 = new QVBoxLayout(layoutWidget);
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setContentsMargins(11, 11, 11, 11);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        verticalLayout_5->setContentsMargins(0, 0, 0, 0);

        label = new QLabel(layoutWidget);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy3(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy3);
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        label->setFont(font);

        label2 = new QLabel(layoutWidget);
        label2->setObjectName(QString::fromUtf8("label2"));
        QSizePolicy sizePolicy6(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label2->setSizePolicy(sizePolicy6);
        QFont font2;
        font2.setBold(true);
        font2.setWeight(75);
        label2->setFont(font2);

        verticalLayout_5->addWidget(label2);

        mpIsLearnCheckBox = new QCheckBox(layoutWidget);
        mpIsLearnCheckBox->setObjectName(QString::fromUtf8("mpIsLearnCheckBox"));

        verticalLayout_5->addWidget(mpIsLearnCheckBox);

        mpNextObjButton = new QPushButton(layoutWidget);
        mpNextObjButton->setObjectName(QString::fromUtf8("mpNextObjButton"));
        QSizePolicy sizePolicy5(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(mpNextObjButton->sizePolicy().hasHeightForWidth());
        mpNextObjButton->setSizePolicy(sizePolicy5);

        verticalLayout_5->addWidget(mpNextObjButton);

        mpPrevObjButton = new QPushButton(layoutWidget);
        mpPrevObjButton->setObjectName(QString::fromUtf8("mpPrevObjButton"));
        QSizePolicy sizePolicy7(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(mpPrevObjButton->sizePolicy().hasHeightForWidth());
        mpPrevObjButton->setSizePolicy(sizePolicy7);

        verticalLayout_5->addWidget(mpPrevObjButton);

       fps_slider = new QSlider( Qt::Horizontal, layoutWidget );
       fps_slider->setObjectName(QString::fromUtf8("mpFpsSlider"));
       fps_slider->setRange( 1, 15 );
       fps_slider->setValue( 7 );
       
       verticalLayout_5->addWidget(fps_slider);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_5->addItem(verticalSpacer);

        mpTabWidget->addTab(Control, QString());
        Flags = new QWidget();
        Flags->setObjectName(QString::fromUtf8("Flags"));
        groupBox_4 = new QGroupBox(Flags);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        groupBox_4->setGeometry(QRect(0, 10, 151, 151));
        layoutWidget_3 = new QWidget(groupBox_4);
        layoutWidget_3->setObjectName(QString::fromUtf8("layoutWidget_3"));
        layoutWidget_3->setGeometry(QRect(1, 20, 169, 108));
        verticalLayout = new QVBoxLayout(layoutWidget_3);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);

        mpShowMessagesCheckBox = new QCheckBox(layoutWidget_3);
        mpShowMessagesCheckBox->setObjectName(QString::fromUtf8("mpShowMessagesCheckBox"));
        mpShowMessagesCheckBox->setChecked(true);

        verticalLayout->addWidget(mpShowMessagesCheckBox);

        mpShowEdgesCheckBox = new QCheckBox(layoutWidget_3);
        mpShowEdgesCheckBox->setObjectName(QString::fromUtf8("mpShowEdgesCheckBox"));

        verticalLayout->addWidget(mpShowEdgesCheckBox);

        mpShowChainsCheckBox = new QCheckBox(layoutWidget_3);
        mpShowChainsCheckBox->setObjectName(QString::fromUtf8("mpShowChainCheckBox"));

        verticalLayout->addWidget(mpShowChainsCheckBox);

        mpShowAllChainsCheckBox = new QCheckBox(layoutWidget_3);
        mpShowAllChainsCheckBox->setObjectName(QString::fromUtf8("mpShowAllChainsCheckBox"));

        verticalLayout->addWidget(mpShowAllChainsCheckBox);

        mpShowLinesCheckBox = new QCheckBox(layoutWidget_3);
        mpShowLinesCheckBox->setObjectName(QString::fromUtf8("mpShowLinesCheckBox"));

        verticalLayout->addWidget(mpShowLinesCheckBox);

        groupBox_5 = new QGroupBox(Flags);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        groupBox_5->setGeometry(QRect(0, 150, 161, 101));
        layoutWidget_4 = new QWidget(groupBox_5);
        layoutWidget_4->setObjectName(QString::fromUtf8("layoutWidget_4"));
        layoutWidget_4->setGeometry(QRect(0, 20, 161, 80));
        verticalLayout_3 = new QVBoxLayout(layoutWidget_4);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);

        mpTabWidget->addTab(Flags, QString());
        MainWindow->setCentralWidget(centralWidget);

        retranslateUi(MainWindow);

        mpTabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Bristol Multi-Obj Detector 2012", 0, QApplication::UnicodeUTF8));
        label2->setText(QApplication::translate("MainWindow", "Learn Objects", 0, QApplication::UnicodeUTF8));
        mpNextObjButton->setText(QApplication::translate("MainWindow", "Next Obj", 0, QApplication::UnicodeUTF8));
        mpPrevObjButton->setText(QApplication::translate("MainWindow", "Prev Obj", 0, QApplication::UnicodeUTF8));
       // fps_slider->setText(QApplication::translate("MainWindow", "FPS", 0, QApplication::UnicodeUTF8));
        mpIsLearnCheckBox->setText(QApplication::translate("MainWindow", "Learning", 0, QApplication::UnicodeUTF8));
        mpTabWidget->setTabText(mpTabWidget->indexOf(Control), QApplication::translate("MainWindow", "Control", 0, QApplication::UnicodeUTF8));
        mpShowMessagesCheckBox->setText(QApplication::translate("MainWindow", "Show Messages (M)", 0, QApplication::UnicodeUTF8));
        mpShowEdgesCheckBox->setText(QApplication::translate("MainWindow", "Show Edges (E)", 0, QApplication::UnicodeUTF8));
        mpShowChainsCheckBox->setText(QApplication::translate("MainWindow", "Show Chains", 0, QApplication::UnicodeUTF8));
        mpShowAllChainsCheckBox->setText(QApplication::translate("MainWindow", "Show All Chains", 0, QApplication::UnicodeUTF8));
        mpShowLinesCheckBox->setText(QApplication::translate("MainWindow", "Show Lines", 0, QApplication::UnicodeUTF8));
        mpTabWidget->setTabText(mpTabWidget->indexOf(Flags), QApplication::translate("MainWindow", "Flags", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
