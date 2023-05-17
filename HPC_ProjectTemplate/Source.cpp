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

void blur(Filter filter, Image_struct input_img, Image_struct* output_img, int offset_rows);
Image_struct add_zero_padding(Image_struct input_img, int filter_size);
void printMatrix(const int* matrix, int m, int n)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << std::setw(5) << matrix[i * n + j] << " ";
		}
		std::cout << std::endl;
	}
}

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
	int* final_data{};

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
	Image_struct padded_img{};
	Image_struct sub_img{};
	Image_struct* final_output_img = new Image_struct();
	int* start_displacement_arr = new int[world_size];
	int* send_counts = new int[world_size];
	int* gather_output_sizes = new int[world_size];
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

		cout << "input Image width = " << input_img.width << "\n" << "input Image height = " << input_img.height << endl;
		
		//add zero padding to input image
		padded_img = add_zero_padding(input_img, filter.size);

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
		int row_displacement = 0;
		start_displacement_arr[0] = row_displacement;
		for (int i = 1; i < world_size; i++) {
			row_displacement = (base_send_count * i) - (offset_rows * sub_img.width);
			start_displacement_arr[i] = row_displacement;
		}
		final_data = new int[input_img.width * input_img.height];
		final_output_img->data = new int[input_img.width * input_img.height];
		final_output_img->width = input_img.width;
		final_output_img->height = input_img.height;
		//free the input image allocated memory
		delete input_img.data;
		//start timer
		start_s = clock();
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
	//broadcast the sendcount array
	MPI_Bcast(send_counts, world_size, MPI_INT, 0, MPI_COMM_WORLD);
	//broadcast the sub image dimensions
	MPI_Bcast(&sub_img.width, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//scatter the padded image
	sub_img.height = send_counts[world_rank] / sub_img.width;
	sub_img.data = new int[((sub_img.height + (filter.size - 1)) * sub_img.width)];//assum the largest size of data
	cout << "Rank "<<world_rank<<" Sub image height = "<<sub_img.height << endl;
	MPI_Scatterv(
		padded_img.data, send_counts, start_displacement_arr, MPI_INT,
		sub_img.data, ((sub_img.height + (filter.size - 1)) * sub_img.width), MPI_INT,
		0, MPI_COMM_WORLD
	);
	//createImage(sub_img.data, sub_img.width, sub_img.height, world_rank);
	

	Image_struct* local_output_img = new Image_struct();

	cout << "Rank  " << world_rank << " sub_img dimensions = " <<sub_img.width<<" x "<<sub_img.height << ", send count = " <<send_counts[world_rank]<< endl;
	blur(filter, sub_img, local_output_img, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	//Gather the local_outputs into the final_output image
	//gather the local_output sizes into gather_displs array
	int local_output_img_size = (local_output_img->width * local_output_img->height);
	int* recieve_displs = new int[world_size];
	MPI_Gather(&local_output_img_size, 1, MPI_INT, gather_output_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (world_rank == 0) {
		//calculate displacements
		recieve_displs[0] = 0;
		for (int i = 1; i < world_size; i++) {
			recieve_displs[i] = 0;
			for (int j = 0; j < i; j++) {
				recieve_displs[i] += gather_output_sizes[j];
			}
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(
		local_output_img->data, (local_output_img->height * local_output_img->width), MPI_INT,
		final_output_img->data, send_counts, recieve_displs,MPI_INT, 0,MPI_COMM_WORLD
		);
		
	if (world_rank == 0) {
		createImage(final_output_img->data, final_output_img->width, final_output_img->height, world_rank + 5);
	
	}
	
	//stop timer
	//stop_s = clock();

	//cout << "output Image width = " << output_img->width << "\n" << "output Image height = " << output_img->height << endl;

	//TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;

	//create the image
	//createImage(output_img->data, output_img->width, output_img->height, world_rank);
	MPI_Barrier(MPI_COMM_WORLD);

	//cout << "time: " << TotalTime << endl;
	delete filter.data;
	delete local_output_img->data;
	delete sub_img.data;
	if (world_rank == 0) {
		delete padded_img.data;
		delete final_output_img->data;
	}
	delete send_counts;
	delete start_displacement_arr;
	delete gather_output_sizes;
	delete recieve_displs;
	MPI_Finalize();

}

void blur(Filter filter, Image_struct input_img, Image_struct* output_img, int offset_rows) {
	output_img->height = (input_img.height - filter.size + 1) - offset_rows;
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
	cout << "Padded image height and width = " << padded_img.height << endl;

	return padded_img;
}