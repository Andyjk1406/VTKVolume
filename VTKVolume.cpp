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
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
 
#include <vtkSphereSource.h>

#include <vtkCommand.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataMapper.h>

#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkMassProperties.h>

// PolyData to process
vtkSmartPointer<vtkPolyData> ClippedpolyData =vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> originalpolyData = vtkSmartPointer<vtkPolyData>::New();

void loadPlanesVec(std::vector<vtkSmartPointer<vtkPlane>>& vPlanes, std::string f_name);
void savePlanesVec(std::vector<vtkSmartPointer<vtkPlane>>& vPlanes, std::string f_name);
void measureVolume(vtkSmartPointer<vtkPolyData> poly);
void saveSTL(vtkSmartPointer<vtkPolyData> poly, std::string fname);
 
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
			else if (key == "c") {

				// Create a sliced model (no volume) from the plane collection
				polyData->DeepCopy(originalpolyData);
				
				//clipper->InsideOutOn();
				
				for (int i = 0; i < vPlanes_.size(); i++) {

					vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
					clipper->SetInputData(polyData);

					clipper->SetClipFunction(vPlanes_[i]);
					clipper->Update();
					ClippedpolyData->DeepCopy(clipper->GetOutput());
					polyData->DeepCopy(ClippedpolyData);
					polyData->Modified();

				}
			
				
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
			else if (key == "1") {

				// Write the vector of planes to a file
				std::string fname;
				std::cout << "Enter filename to write planes to : ";
				std::cin >> fname;
				savePlanesVec(vPlanes_, fname);

			}
			else if (key == "2") {

				// Read the vector of planes to a file
				std::string fname;
				std::cout << "Enter filename to read planes from : ";
				std::cin >> fname;
				loadPlanesVec(vPlanes_, fname);

			}
			else if (key == "v") {

				measureVolume(polyData);

			}
			else if (key == "3") {

				// Save the STL
				std::string fname;
				std::cout << "Enter filename to save PLY (without the '.ply' extension) : ";
				std::cin >> fname;
				saveSTL(polyData, fname);

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

	if (argc == 5) // .exe file_input plane_crop action output 
	{
		// Batch script run

		// Load the mesh
		vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		polyData->DeepCopy(reader->GetOutput());
		originalpolyData->DeepCopy(reader->GetOutput());

		// Load the planes file
		std::vector<vtkSmartPointer<vtkPlane>> vPlanes;
		std::string planes_file = argv[2];
		planes_file += ".pln";
		loadPlanesVec(vPlanes, planes_file);

		// Do the cut/fill
		if (argv[3] != "fill") // Any text other than 'fill' will induce a cut
		{
			for (int i = 0; i < vPlanes.size(); i++) {

				vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
				clipper->SetInputData(polyData);
				clipper->SetClipFunction(vPlanes[i]);
				clipper->Update();
				polyData->DeepCopy(clipper->GetOutput());
				polyData->Modified();

			}
		}
		else //fill
		{
			// Create the volume polydata
			vtkSmartPointer<vtkPlaneCollection> planes = vtkSmartPointer<vtkPlaneCollection>::New();
			for (int i = 0; i<vPlanes.size(); i++)
				planes->AddItem(vPlanes[i]);


			vtkSmartPointer<vtkClipClosedSurface> clipper = vtkSmartPointer<vtkClipClosedSurface>::New();
			clipper->SetInputData(originalpolyData);
			clipper->SetClippingPlanes(planes);
			clipper->SetActivePlaneId(2);
			clipper->SetScalarModeToColors();
			clipper->SetClipColor(0.8900, 0.8100, 0.3400); // banana
			clipper->SetBaseColor(1.0000, 0.3882, 0.2784); // tomato
			clipper->SetActivePlaneColor(0.6400, 0.5800, 0.5000); // beige
			clipper->Update();

			polyData->DeepCopy(clipper->GetOutput());
			polyData->Modified();

		}

		// Save the resulting file 
		std::string save_name = argv[4];
		saveSTL(polyData, save_name); // NOTE! This actually saves PLYs at the moment (including tacking on the .ply)!

		// Quit
		return 1;

	}

  else if (argc > 1)
    {
    vtkSmartPointer<vtkPLYReader> reader = 
      vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    polyData->DeepCopy(reader->GetOutput());
	originalpolyData->DeepCopy(reader->GetOutput());
    }
  else
    {

	  cout << "Takes one PLY file" << endl;
	  cout << "s - add slice, f - create filled block, 1 - save slicing planes, 2 - load slicing planes, 3 - Save the STL" << endl;
	  cout << "v - measure volume, c - clip with a loaded set of planes but no volume filling" << endl;
	  cout << "or '.exe file_input plane_crop action(fill/cut) output'" << endl;

	  // Create a sphere
    vtkSmartPointer<vtkSphereSource> sphereSource = 
      vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetThetaResolution(20);
    sphereSource->SetPhiResolution(11);
    sphereSource->Update();
 
    polyData = sphereSource->GetOutput();
    }


  // Setup a visualization pipeline
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
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
  
  return EXIT_SUCCESS;
}

void savePlanesVec(std::vector<vtkSmartPointer<vtkPlane>>& vPlanes, std::string f_name) {

	ofstream myFile(f_name+".pln", ios::binary);

	if (myFile.is_open()) {
		int n = vPlanes.size();
		myFile.write(reinterpret_cast<const char *>(&n), sizeof(n));

		for (int i = 0; i < n; i++) {
			double norm[3];
			double orig[3];
			vPlanes[i]->GetNormal(norm);
			vPlanes[i]->GetOrigin(orig);
			myFile.write(reinterpret_cast<const char *>(&orig), sizeof(orig));
			myFile.write(reinterpret_cast<const char *>(&norm), sizeof(norm));
		}
	}
	myFile.close();

	cout << "File saved..." << endl;

	for (int i = 0; i < vPlanes.size(); i++) {
		double norm[3];
		double orig[3];
		vPlanes[i]->GetNormal(norm);
		vPlanes[i]->GetOrigin(orig);
		cout << "Origin : " << orig[0] << "," << orig[1] << "," << orig[2] << endl;
		cout << "Normal : " << norm[0] << "," << norm[1] << "," << norm[2] << endl;
	}

}


void loadPlanesVec(std::vector<vtkSmartPointer<vtkPlane>>& vPlanes, std::string f_name) {

	vPlanes.clear();
	ifstream myFile(f_name+".pln", ios::binary);

	if (myFile.is_open()) {
		int n=0;
		myFile.read(reinterpret_cast<char *>(&n), sizeof(n));

		for (int i = 0; i < n; i++) {
			double norm[3];
			double orig[3];
			myFile.read(reinterpret_cast<char *>(&orig), sizeof(orig));
			myFile.read(reinterpret_cast<char *>(&norm), sizeof(norm));
			vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
			plane->SetOrigin(orig);
			plane->SetNormal(norm);
			vPlanes.push_back(plane);
		}

	}
	myFile.close();

	cout << "File loaded..." << endl;

	for (int i = 0; i < vPlanes.size(); i++) {
		double norm[3];
		double orig[3];
		vPlanes[i]->GetNormal(norm);
		vPlanes[i]->GetOrigin(orig);
		cout << "Origin : " << orig[0] << "," << orig[1] << "," << orig[2] << endl;
		cout << "Normal : " << norm[0] << "," << norm[1] << "," << norm[2] << endl;
	}
}

void measureVolume(vtkSmartPointer<vtkPolyData> poly) {

	vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPolyData->SetInputData(poly);
	cleanPolyData->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputConnection(cleanPolyData->GetOutputPort());
	triangleFilter->Update();

	vtkSmartPointer<vtkMassProperties> mass = vtkSmartPointer<vtkMassProperties>::New();
	mass->SetInputConnection(triangleFilter->GetOutputPort());
	mass->Update();

	double volume = mass->GetVolume();
	cout << "Mesh volume : " << volume << endl;

}

void saveSTL(vtkSmartPointer<vtkPolyData> poly, std::string fname)
{
	// Save
	vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileTypeToBinary();
	std::string STLfname = fname + ".ply";
	writer->SetFileName(STLfname.c_str());
	writer->Write();
}