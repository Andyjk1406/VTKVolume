#include <vtkVersion.h>
#include <vtkSmartPointer.h>
 
#include <vtkActor.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPLYReader.h>
#include <vtkCleanPolyData.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>
#include <vtkSphereSource.h>
#include <vtkBooleanOperationPolyDataFilter.h>

float standard_deviation(vtkSmartPointer<vtkDistancePolyDataFilter> dist);
float test;
float finalTestGIT;

//Help!
//Test comment and a new bit

std::string testFilename;
 
int main(int argc, char* argv[])
{
  vtkSmartPointer<vtkPolyData> input1;
  vtkSmartPointer<vtkPolyData> input2;
  if (argc == 3)
    {
//     std::cerr << "Usage: " << argv[0]
//               << " filename1.ply"
//               << " filename2.ply" << std::endl;
		testFilename = argv[1];

vtkSmartPointer<vtkPLYReader> reader1 = vtkSmartPointer<vtkPLYReader>::New();
	reader1->SetFileName(argv[1]);
	reader1->Update();
	input1 = reader1->GetOutput();	

vtkSmartPointer<vtkPLYReader> reader2 = vtkSmartPointer<vtkPLYReader>::New();
	reader2->SetFileName(argv[2]);
	reader2->Update();
	input2 = reader2->GetOutput();		
	
    
    }
  else
    {
    vtkSmartPointer<vtkSphereSource> sphereSource1 =
      vtkSmartPointer<vtkSphereSource>::New();
    sphereSource1->SetCenter(1, 0, 0);
    sphereSource1->Update();
    input1 = sphereSource1->GetOutput();
 
    vtkSmartPointer<vtkSphereSource> sphereSource2 =
      vtkSmartPointer<vtkSphereSource>::New();
    sphereSource2->Update();
    input2 = sphereSource2->GetOutput();
    }
 
  vtkSmartPointer<vtkCleanPolyData> clean1 =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  clean1->SetInputConnection( input1->GetProducerPort());
#else
  clean1->SetInputData( input1);
#endif
 
  vtkSmartPointer<vtkCleanPolyData> clean2 =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  clean2->SetInputConnection( input2->GetProducerPort());
#else
  clean2->SetInputData( input2);
#endif
 
  vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter =
    vtkSmartPointer<vtkDistancePolyDataFilter>::New();
 
  distanceFilter->SetInputConnection( 0, clean1->GetOutputPort() );
  distanceFilter->SetInputConnection( 1, clean2->GetOutputPort() );
  distanceFilter->Update();
/*
  // Boolean
  vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOp = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
  booleanOp->SetOperationToIntersection();
  booleanOp->SetInputConnection(0, clean1->GetOutputPort());
  booleanOp->SetInputConnection(1, clean2->GetOutputPort());

  vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  bmapper->SetInputConnection(booleanOp->GetOutputPort());
  bmapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> bactor =vtkSmartPointer<vtkActor>::New();
  bactor->SetMapper(bmapper);
*/

  float sd = standard_deviation(distanceFilter);
 
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( distanceFilter->GetOutputPort() );
  mapper->SetScalarRange(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetRange()[0],
   distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetRange()[1]);
 
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper( mapper );
 
  vtkSmartPointer<vtkScalarBarActor> scalarBar = 
    vtkSmartPointer<vtkScalarBarActor>::New();
  scalarBar->SetLookupTable(mapper->GetLookupTable());
  scalarBar->SetTitle("Distance");
  scalarBar->SetNumberOfLabels(4);
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
 
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer( renderer );
 
  vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renWinInteractor->SetRenderWindow( renWin );
 
  renderer->AddActor( actor );
 // renderer->AddActor(bactor);
  renderer->AddActor2D(scalarBar);
 
  //renWin->Render();
  //renWinInteractor->Start();
 
  return EXIT_SUCCESS;
}

float standard_deviation(vtkSmartPointer<vtkDistancePolyDataFilter> dist)
{
	double meanNeg = 0.0, sum_deviationNeg = 0.0;
	double meanPos = 0.0, sum_deviationPos = 0.0;
	double mean = 0.0, sum_deviation = 0.0;
	double unsigned_mean = 0.0, unsigned_sum_dev = 0.0;
	vtkIdType nValues;
	int meanNegctr = 0, meanPosctr = 0;
	nValues = dist->GetOutput()->GetPointData()->GetScalars()->GetSize();

	for (vtkIdType i = 0; i<nValues; ++i)
	{
		double val = dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
		mean += val;
		if (val < 0) {	meanNeg += val; meanNegctr++;}
		if (val >= 0) { meanPos += val; meanPosctr++; }
		unsigned_mean += fabs(dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i));
	}

	cout << "MeanNeg : " << meanNeg << endl;
	cout << "Neg count : " << meanNegctr << endl;

	mean = mean / (double)nValues;
	meanNeg = meanNeg / (double)meanNegctr;
	meanPos = meanPos / (double)meanPosctr;
	unsigned_mean = unsigned_mean / (double)nValues;

	double negCtr = 0; double posCtr = 0;

	for (vtkIdType i = 0; i < nValues; ++i) {
		sum_deviation += (dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - mean)*(dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - mean);
		if (dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) < 0) {
			sum_deviationNeg += (dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - meanNeg)*(dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - meanNeg);
			negCtr++;
		}
		else {
			sum_deviationPos += (dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - meanPos)*(dist->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i) - meanPos);
			posCtr++;
		}
	}
		
	sum_deviation = sqrt(sum_deviation / (double)nValues);
	sum_deviationNeg = sqrt(sum_deviationNeg / negCtr);
	sum_deviationPos = sqrt(sum_deviationPos / posCtr);

	ofstream out("DataResults.txt", ios::app);
	out << testFilename;
	out << " Mean " << mean;
	out << " UMean " << unsigned_mean;
	out << " Mean_Neg " << meanNeg;
	out << " Mean_Pos " << meanPos;
	out << " SD " << sum_deviation;
	out << " SDpos " << sum_deviationPos;
	out << " SDneg " << sum_deviationNeg;
	out << " Range " << dist->GetOutput()->GetPointData()->GetScalars()->GetRange()[0] << "  " << dist->GetOutput()->GetPointData()->GetScalars()->GetRange()[1] << std::endl;

	out.close();


	std::cout << "Mean " << mean;
	std::cout << " UMean " << unsigned_mean;
	std::cout << " Mean_Neg " << meanNeg;
	std::cout << " Mean_Pos " << meanPos;
	std::cout << " SD " << sum_deviation;
	std::cout << " SDpos " << sum_deviationPos;
	std::cout << " SDneg " << sum_deviationNeg;
	std::cout << " Range " << dist->GetOutput()->GetPointData()->GetScalars()->GetRange()[0] << "  " << dist->GetOutput()->GetPointData()->GetScalars()->GetRange()[1] << std::endl;


	return sum_deviation;
}