#include<mpi.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>

using namespace std;
using namespace msclr::interop;

struct Filter
{
	int size;
	int* data;
};
struct Image_struct {
	int width;
	int height;
	int* data;
};
#include <iostream>
#include <iomanip>

void blur(Filter filter, Image_struct input_img, Image_struct* output_img);
Image_struct add_zero_padding(Image_struct input_img, int filter_size);

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int* Red = new int[BM.Height * BM.Width];
	int* Green = new int[BM.Height * BM.Width];
	int* Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height * BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i * BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	cout << "result Image Saved " << index << endl;
}

int main()
{

	MPI_Init(NULL, NULL);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int start_s, stop_s, TotalTime = 0;

	//Image Path
	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//lena.png";

	Filter filter;
	Image_struct padded_img{};	// will be scattered among processes
	Image_struct sub_img{};		// will hold the scattered data
	Image_struct* final_output_img = new Image_struct();	// will hold the gathered results from processes
	int* start_displacement_arr{};							// Entry i specifies the displacement (relative to the send-buffer) from which to take the outgoing padded image data to process i
	int* send_counts = new int[world_size];					// Entry i specifies the number of elements to send to each process from the send-buffer
	int* gather_output_sizes{};								// Holds the output size of each process so that Gatherv can know how much to receive from each process
	int* recieve_displs{};									// Entry i specifies the displacement relative to recvbuf at which to place the incoming data from process i into the final output image
	if (world_rank == 0) {
		//Filter intialization
		cout << "Enter kernel size (must be an odd number)" << endl;
		cin >> filter.size;
		if (filter.size % 2 == 0)
			return 0;
		filter.data = new int[filter.size * filter.size];
		for (int i = 0; i < filter.size; i++) {
			for (int j = 0; j < filter.size; j++) {
				filter.data[i * filter.size + j] = 1;
			}
		}
		//input the image
		Image_struct input_img;
		imagePath = marshal_as<System::String^>(img);
		input_img.data = inputImage(&input_img.width, &input_img.height, imagePath);

		cout << "input Image dimensions = " << input_img.width << " x " << input_img.height << endl;

		//add zero padding to input image
		padded_img = add_zero_padding(input_img, filter.size);
		cout << "padded Image dimensions = " << padded_img.width << " x " << padded_img.height << endl;

		//the timer starts exactly before the code written for MPI runs
		start_s = clock();
		
		// the image will be split at the rows with full width
		sub_img.height = (padded_img.height / world_size);
		sub_img.width = padded_img.width; 

		//calculate how much to send from the padded image to each process (all of them get same count except for 0)
		int base_send_count = sub_img.height * sub_img.width;
		int offset_rows = filter.size - 1;
		send_counts[0] = base_send_count;
		for (int i = 1; i < world_size; i++) {
			send_counts[i] = base_send_count + (offset_rows * sub_img.width);
		}
		//calculate the start of every process in the padded image
		start_displacement_arr = new int[world_size];
		int displacement = 0;
		start_displacement_arr[0] = displacement;
		for (int i = 1; i < world_size; i++) {
			displacement = (base_send_count * i) - (offset_rows * sub_img.width);
			start_displacement_arr[i] = displacement;
		}
		//setup the final output image
		final_output_img->data = new int[input_img.width * input_img.height];
		final_output_img->width = input_img.width;
		final_output_img->height = input_img.height;

		//We won't need the input image again in the code
		delete input_img.data;

		// these are only needed by procees 0
		gather_output_sizes = new int[world_size];
		recieve_displs = new int[world_size];

		//broadcast the filter
		MPI_Bcast(&filter.size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(filter.data, (filter.size * filter.size), MPI_INT, 0, MPI_COMM_WORLD);

	}
	else {
		//recieve the filter ?
		MPI_Bcast(&filter.size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		filter.data = new int[filter.size * filter.size];
		MPI_Bcast(filter.data, (filter.size * filter.size), MPI_INT, 0, MPI_COMM_WORLD);

	}
	// Each process uses the send_counts array to calculate the height of it's sub image
	MPI_Bcast(send_counts, world_size, MPI_INT, 0, MPI_COMM_WORLD);
	// All processes will have the same sub image width
	MPI_Bcast(&sub_img.width, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// calculate the height based on how much the process should receive
	sub_img.height = send_counts[world_rank] / sub_img.width;
	sub_img.data = new int[send_counts[world_rank]];
	
	//scatter the padded image
	MPI_Scatterv(
		padded_img.data, send_counts, start_displacement_arr, MPI_INT,
		sub_img.data, (sub_img.height * sub_img.width), MPI_INT,
		0, MPI_COMM_WORLD
	);	

	Image_struct* local_output_img = new Image_struct();

	cout << "Rank " << world_rank << " sub image dimensions = " <<sub_img.width<<" x "<<sub_img.height << ", send count = " <<send_counts[world_rank]<< endl;
	blur(filter, sub_img, local_output_img);

	//Gathers all the local-output-images sizes so that they can be used to calculate the displacement of each local-output-image into the final image
	int local_output_img_size = (local_output_img->width * local_output_img->height);
	MPI_Gather(&local_output_img_size, 1, MPI_INT, gather_output_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (world_rank == 0) {
		// each local-output-image is placed at a displacement that is equal to the sum of the sizes of the local-output-images before it
		recieve_displs[0] = 0;
		for (int i = 1; i < world_size; i++) {
			recieve_displs[i] = 0;
			for (int j = 0; j < i; j++) {
				recieve_displs[i] += gather_output_sizes[j];
			}
		}

	}

	// Gather the final-output-image from all the local-output-images
	MPI_Gatherv(
		local_output_img->data, (local_output_img->height * local_output_img->width), MPI_INT,
		final_output_img->data, send_counts, recieve_displs,MPI_INT, 0,MPI_COMM_WORLD
		);
	// process 0 is responsible for creating the final image
	if (world_rank == 0) {
		// the timer stops right after the final output image is gathered from all the processes
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		cout << "time: " << TotalTime << endl;
		
		cout << "output Image dimensions = " << final_output_img->width << " x " << final_output_img->height << endl;
		createImage(final_output_img->data, final_output_img->width, final_output_img->height, world_rank);
	}
	
	delete filter.data;
	delete local_output_img->data;
	delete sub_img.data;
	delete send_counts;
	if (world_rank == 0) {
		delete padded_img.data;
		delete final_output_img->data;
		delete start_displacement_arr;
		delete gather_output_sizes;
		delete recieve_displs;

	}
	MPI_Finalize();

}

void blur(Filter filter, Image_struct input_img, Image_struct* output_img) {
	output_img->height = (input_img.height - filter.size + 1);
	output_img->width = input_img.width - filter.size + 1;

	output_img->data = new int[output_img->height * output_img->width];
	int sop = 0;
	for (int i = 0; i < output_img->height ; i++) {
		for (int j = 0; j < output_img->width; j++) {
			sop = 0;
			for (int h = 0; h < filter.size; h++) {
				for (int w = 0; w < filter.size; w++) {
					sop += filter.data[h * filter.size + w] * input_img.data[(i + h) * input_img.width + (j + w)];
				}
			}
			output_img->data[i *output_img->width + j] = sop / (filter.size * filter.size);
		}
	}

}
Image_struct add_zero_padding(Image_struct input_img, int filter_size) {
	Image_struct padded_img;
	padded_img.height = input_img.height + filter_size - 1;
	padded_img.width = input_img.width + filter_size - 1;
	padded_img.data = new int[padded_img.height * padded_img.width];
	/* by default the padded image is initialized with zeros. We just have to fill in with the input image */
	int input_img_start_row = (filter_size - 1) / 2;//where the input image should start inside the padded image
	int input_img_start_col = (filter_size - 1) / 2;
	for (int i = 0; i < input_img.height; i++) {
		for (int j = 0; j < input_img.width; j++) {
			padded_img.data[(i + input_img_start_row) * padded_img.width + (j + input_img_start_col)]
				= input_img.data[i * input_img.width + j];
		}
	}
	return padded_img;
}