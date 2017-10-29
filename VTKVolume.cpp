#include <vtkVersion.h>
#include <vtkSmartPointer.h>
 
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkClipClosedSurface.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
 
#include <vtkSphereSource.h>

#include <vtkCommand.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataMapper.h>

// PolyData to process
vtkSmartPointer<vtkPolyData> ClippedpolyData =vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> originalpolyData = vtkSmartPointer<vtkPolyData>::New();

 
//
// Demonstrate the use of clipping of polygonal data and use of implicitplanewidget2
//

// Callback for the interaction
// This does the actual work: updates the vtkPlane implicit function.
// This in turn causes the pipeline to update and clip the object.
class vtkIPWCallback : public vtkCommand
{
public:
	static vtkIPWCallback *New()
	{
		return new vtkIPWCallback;
	}
	virtual void Execute(vtkObject *caller, unsigned long event, void* calldata)
	{
		if (event == vtkCommand::KeyPressEvent)
		{
			//std::cout << "key press" << std::endl;
			// Get the keypress
			vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
			std::string key = iren->GetKeySym();

			// Output the key that was pressed
			//std::cout << "Pressed " << key << std::endl;

			if (key == "s") {

				// Add the plane to the vector and save the currently sliced polydata			
				clipPoly->Update();
				ClippedpolyData->DeepCopy(clipPoly->GetOutput());
				polyData->DeepCopy(ClippedpolyData);
				polyData->Modified();

				vtkSmartPointer<vtkPlane> currentPlane = vtkSmartPointer<vtkPlane>::New();
				double norm[3];
				double orig[3];
				Plane->GetNormal(norm);
				for (int i = 0; i < 3; i++)
					norm[i] *= -1.;
				Plane->GetOrigin(orig);
				currentPlane->SetNormal(norm);
				currentPlane->SetOrigin(orig);
				vPlanes_.push_back(currentPlane);

			}
			else if (key == "f") {

				// Create the volume polydata
				vtkSmartPointer<vtkPlaneCollection> planes = vtkSmartPointer<vtkPlaneCollection>::New();
				for (int i=0;i<vPlanes_.size();i++)
					planes->AddItem(vPlanes_[i]);
				

				vtkSmartPointer<vtkClipClosedSurface> clipper = vtkSmartPointer<vtkClipClosedSurface>::New();
#if VTK_MAJOR_VERSION <= 5
				clipper->SetInput(originalpolyData);
#else
				clipper->SetInputData(originalpolyData);
#endif
				clipper->SetClippingPlanes(planes);
				clipper->SetActivePlaneId(2);
				clipper->SetScalarModeToColors();
				clipper->SetClipColor(0.8900, 0.8100, 0.3400); // banana
				clipper->SetBaseColor(1.0000, 0.3882, 0.2784); // tomato
				clipper->SetActivePlaneColor(0.6400, 0.5800, 0.5000); // beige
				clipper->Update();

				ClippedpolyData->DeepCopy(clipper->GetOutput());
				polyData->DeepCopy(ClippedpolyData);
				polyData->Modified();

			}
			
		}
		else {


			vtkImplicitPlaneWidget2 *planeWidget =
				reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
			vtkImplicitPlaneRepresentation *rep =
				reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
			rep->GetPlane(this->Plane);
		}

	}
	vtkIPWCallback() :Plane(0), Actor(0), clipPoly(0) {}
	vtkPlane *Plane;
	vtkActor *Actor;
	vtkClipPolyData *clipPoly;
	std::vector<vtkSmartPointer<vtkPlane>> vPlanes_;

};
 
int main (int argc, char *argv[])
{
	//vtkSmartPointer<vtkPolyData> polyData;
  if (argc > 1)
    {
    vtkSmartPointer<vtkSTLReader> reader = 
      vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    polyData->DeepCopy(reader->GetOutput());
	originalpolyData->DeepCopy(reader->GetOutput());
    }
  else
    {
    // Create a sphere
    vtkSmartPointer<vtkSphereSource> sphereSource = 
      vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetThetaResolution(20);
    sphereSource->SetPhiResolution(11);
    sphereSource->Update();
 
    polyData = sphereSource->GetOutput();
    }


  // Setup a visualization pipeline
  vtkSmartPointer<vtkPlane> plane =
	  vtkSmartPointer<vtkPlane>::New();
  vtkSmartPointer<vtkClipPolyData> clipper2 = vtkSmartPointer<vtkClipPolyData>::New();

  
  clipper2->SetClipFunction(plane);
  clipper2->InsideOutOn();
  clipper2->SetInputData(polyData);


  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper =
	  vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(clipper2->GetOutputPort());
  vtkSmartPointer<vtkActor> actor =
	  vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkProperty> backFaces =
	  vtkSmartPointer<vtkProperty>::New();
  backFaces->SetDiffuseColor(.8, .8, .4);

  actor->SetBackfaceProperty(backFaces);

  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
	  vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
	  vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderer->AddActor(actor);

  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderWindow->Render();

  // The callback will do the work
  vtkSmartPointer<vtkIPWCallback> myCallback =
	  vtkSmartPointer<vtkIPWCallback>::New();
  myCallback->Plane = plane;
  myCallback->Actor = actor;
  myCallback->clipPoly = clipper2;

  vtkSmartPointer<vtkImplicitPlaneRepresentation> rep =
	  vtkSmartPointer<vtkImplicitPlaneRepresentation>::New();
  rep->SetPlaceFactor(1.25); // This must be set prior to placing the widget
  rep->PlaceWidget(actor->GetBounds());
  rep->SetNormal(plane->GetNormal());

  vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget =
	  vtkSmartPointer<vtkImplicitPlaneWidget2>::New();
  planeWidget->SetInteractor(renderWindowInteractor);
  planeWidget->SetRepresentation(rep);
  planeWidget->AddObserver(vtkCommand::InteractionEvent, myCallback);
  renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, myCallback);

  // Render

  renderWindowInteractor->Initialize();
  renderWindow->Render();
  planeWidget->On();

  // Begin mouse interaction
  renderWindowInteractor->Start();
  
 /*
  double *center = polyData->GetCenter();
  vtkSmartPointer<vtkPlane> plane1 =
    vtkSmartPointer<vtkPlane>::New();
  plane1->SetOrigin(center[0], center[1], center[2]);
  plane1->SetNormal(0.0, -1.0, 0.0);
  vtkSmartPointer<vtkPlane> plane2 =
    vtkSmartPointer<vtkPlane>::New();
  plane2->SetOrigin(center[0], center[1], center[2]);
  plane2->SetNormal(0.0, 0.0, 1.0);
  vtkSmartPointer<vtkPlane> plane3 =
    vtkSmartPointer<vtkPlane>::New();
  plane3->SetOrigin(center[0], center[1], center[2]);
  plane3->SetNormal(-1.0, 0.0, 0.0);

  // Cecilies custom planes
  //6.9893404x - 1.925634y + 0.21222094z - 91.87394653458 = 0
  //PP0; 45.1051; 94.6566; -193.699
  plane1->SetOrigin(45.1051, 94.6566, -193.699);
  plane1->SetNormal(-6.9893404, 1.925634, -0.21222094);
  // - 0.3512672x + 5.8355029y - 0.17647366z - 581.38209281627 = 0
  //PP0; 48.675; 96.6949; -193.891
  plane2->SetOrigin(48.675, 96.6949, -193.891);
  plane2->SetNormal(-0.3512672, 5.8355029,  -0.17647366);
  // 7.9424392x + 3.4904806y - 0.26454562z - 797.61333832214 = 0
  // PP0; 52.3309; 94.6747; -194.742
  plane3->SetOrigin(52.3309, 94.6747, -194.742);
  plane3->SetNormal(-7.9424392, -3.4904806, 0.26454562);
  // 4.4600957x - 2.2531021y - 0.15174307z - 57.01920096349 = 0
  // PP0; 51.0764; 88.9237; -194.854
  vtkSmartPointer<vtkPlane> plane4 =
	  vtkSmartPointer<vtkPlane>::New();
  plane4->SetOrigin(51.0764, 88.9237, -194.854);
  plane4->SetNormal(-4.4600957,  2.2531021,  0.15174307);
  // 3.113502x + 6.23406y + 0.2050836z - 654.811598694 = 0
  // PP0; 45.9838; 88.4242; -193.096
  vtkSmartPointer<vtkPlane> plane5 =
	  vtkSmartPointer<vtkPlane>::New();
  plane5->SetOrigin(45.9838, 88.4242, -193.096);
  plane5->SetNormal(3.113502, 6.23406, 0.2050836);

 
  vtkSmartPointer<vtkPlaneCollection> planes =
    vtkSmartPointer<vtkPlaneCollection>::New();
  planes->AddItem(plane1);
  planes->AddItem(plane2);
  planes->AddItem(plane3);
  planes->AddItem(plane4);
  planes->AddItem(plane5);
 
  vtkSmartPointer<vtkClipClosedSurface> clipper =
    vtkSmartPointer<vtkClipClosedSurface>::New();
#if VTK_MAJOR_VERSION <= 5
  clipper->SetInput(polyData);
#else
  clipper->SetInputData(polyData);
#endif
  clipper->SetClippingPlanes(planes);
  clipper->SetActivePlaneId(2);
  clipper->SetScalarModeToColors();
  clipper->SetClipColor(0.8900, 0.8100, 0.3400); // banana
  clipper->SetBaseColor(1.0000, 0.3882, 0.2784); // tomato
  clipper->SetActivePlaneColor(0.6400, 0.5800, 0.5000); // beige
 
  vtkSmartPointer<vtkDataSetMapper> clipMapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  clipMapper->SetInputConnection(clipper->GetOutputPort());
 
  vtkSmartPointer<vtkActor> clipActor =
    vtkSmartPointer<vtkActor>::New();
  clipActor->SetMapper(clipMapper);
  clipActor->GetProperty()->SetColor(1.0000,0.3882,0.2784);
  clipActor->GetProperty()->SetInterpolationToFlat();
 
  // Create graphics stuff
  //
  vtkSmartPointer<vtkRenderer> ren1 =
    vtkSmartPointer<vtkRenderer>::New();
  ren1->SetBackground(.3, .4, .6);
 
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  renWin->SetSize(512,512);
 
  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);
 
  // Add the actors to the renderer, set the background and size
  //
  ren1->AddActor(clipActor);
 
  // Generate an interesting view
  //
  ren1->ResetCamera();
  ren1->GetActiveCamera()->Azimuth(120);
  ren1->GetActiveCamera()->Elevation(30);
  ren1->GetActiveCamera()->Dolly(1.0);
  ren1->ResetCameraClippingRange();
 
  iren->Initialize();
  iren->Start();
 */
  return EXIT_SUCCESS;
}
