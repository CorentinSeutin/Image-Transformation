package pdl.backend.AlgorithmImage;

import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;

import boofcv.struct.border.BorderType;
import boofcv.io.image.UtilImageIO;
import boofcv.struct.image.GrayU8;
import boofcv.struct.image.Planar;
import boofcv.alg.color.ColorHsv;
import boofcv.alg.filter.binary.GThresholdImageOps;
import boofcv.alg.filter.binary.ThresholdImageOps;
import boofcv.io.image.ConvertBufferedImage;

public class Color {
////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////TP3////////////////////////////////////////////////////////
    public static void brightnessImage(Planar<GrayU8> input, int delta) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

		for (int i = 0; i < input.getNumBands(); ++i){
            final int i_bis = i;

            executor.submit(() -> {
                GrayU8 inputBand = input.getBand(i_bis);
                GrayLevelProcessing.brightness(inputBand, delta);
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	}
    public static void brightnessImageParallel(Planar<GrayU8> input, int delta) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        List<Callable<Void>> tasks = new ArrayList<>();
        for (int i = 0; i < input.getNumBands(); ++i) {
            final int i_bis = i;
            tasks.add(() -> {
                GrayU8 inputBand = input.getBand(i_bis);
                GrayLevelProcessing.brightness(inputBand, delta);
                return null;
            });
        }

        executor.invokeAll(tasks);
        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

    }
////////////////////////////////////////////////////////////////////////////////////////////

	public static int min(GrayU8 input){
		int min = 255;
		for (int y = 0; y < input.height; ++y) {
			for (int x = 0; x < input.width; ++x) {
				int gl = input.get(x, y);
				if (min > gl){
					min = gl;
				}
			}
		}
		return min;
	}

	public static int max(GrayU8 input){
		int max = 0;
		for (int y = 0; y < input.height; ++y) {
			for (int x = 0; x < input.width; ++x) {
				int gl = input.get(x, y);
				if(max < gl){
					max = gl;
				}
			}
		}
		return max;
	}

	public static int lut(GrayU8 input, int min, int max, int ng){
		return (int)( (255*(ng-min))/(max-min) );
	}

	public static void transformLineaire(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int i = 0; i < image.getNumBands(); ++i){

            final int i_bis = i;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(i_bis);
                final GrayU8 output = imgout.getBand(i_bis);

                int min = min(input);
                int max = max(input);

                for (int y = 0; y < input.height; ++y) {
                    for (int x = 0; x < input.width; ++x) {
                        if(i_bis == 3){
                            output.set(x,y,input.get(x,y));
                            continue;
                        }

                        int gl = input.get(x, y);
                        output.set(x, y, lut(input, min, max, gl));
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	}

////////////////////////////////////////////////////////////////////////////////////////////

    public static int lutGeneral(GrayU8 input, int min1, int max1, int min2, int max2, int ng){
        return (int)( ( (ng - min1) * (max2 - min2)/(max1 - min1) ) + min2 );
    }

    public static void transformLineaireMinMax(Planar<GrayU8> image, Planar<GrayU8> imgout, int min2, int max2) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int i = 0; i < image.getNumBands(); ++i){

            final int i_bis = i;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(i_bis);
                final GrayU8 output = imgout.getBand(i_bis);

                int min1 = min(input);
                int max1 = max(input);

                for (int y = 0; y < input.height; ++y) {
                    for (int x = 0; x < input.width; ++x) {
                        if(i_bis == 3){
                            output.set(x,y,input.get(x,y));
                            continue;
                        }

                        int gl = input.get(x, y);
                        output.set(x, y, lutGeneral(input, min1, max1, min2, max2, gl));
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

////////////////////////////////////////////////////////////////////////////////////////////

	public static int[] lutTable(GrayU8 input, int min, int max) {
		int[] tab = new int[256];

		for(int cpt = 0; cpt < 256; cpt++){
            final int cpt_bis = cpt;
			tab[cpt_bis] = lut(input, min, max,cpt_bis);
		}

		return tab;
	}

	public static void transformLineaireLUTtab(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int i = 0; i < image.getNumBands(); ++i){
            final int i_bis = i;

            executor.submit(() -> {
                final GrayU8 input = image.getBand(i_bis);
                final GrayU8 output = imgout.getBand(i_bis);

                int min = min(input);
                int max = max(input);
                int lutTab[];
                lutTab = lutTable(input, min, max);

                for (int y = 0; y < input.height; ++y) {
                    for (int x = 0; x < input.width; ++x) {
                        if(i_bis == 3){
                            output.set(x,y,input.get(x,y));
                            continue;
                        }

                        int gl = input.get(x, y);
                        output.set(x, y, lutTab[gl]);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////TP4/////////////////////////////////////////////////////

    public static int[][] updateMat(GrayU8 input, int i, int j, int n){
        //y ligne x colonne
        int m[][] = new int[(2*n)+1][(2*n)+1];

        for(int y = i-n; y <= i+n; y++){
            final int y_bis = y;
            for(int x = j-n; x <= j+n; x++){
                m[y_bis - (i-n)][x - (j-n)] = input.get(x, y_bis);
            }
        }
        return m;
    }

    public static void updatePixel(GrayU8 input, int m[][], int x, int y, int n){
        int sum = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                sum += m[i][j];
            }
        }

        int nb_pixel = ((2*n) + 1) * ((2*n) + 1);
        sum /= nb_pixel;
        input.set(x, y, sum);
    }

    public static void meanFilterSimple(Planar<GrayU8> image, Planar<GrayU8> imgout, int size) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int n = (size - 1)/2;

                for(int i = n; i < input.height - n; i++) {
                    for(int j = n; j < input.width - n; j++) {
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        int m[][];
                        m = updateMat(input,i,j,n);
                        updatePixel(output,m,j,i,n);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static int[][] updateMatNormalized(GrayU8 input, int i, int j, int n) {
        //y ligne x colonne
        int m[][] = new int[(2*n)+1][(2*n)+1];
        for(int y = i-n; y <= i+n; y++){
            final int y_bis = y;
            for(int x = j-n; x <= j+n; x++){
                if(x < 0 || y_bis < 0 || y_bis >= input.height || x >= input.width){
                    m[y_bis - (i-n)][x - (j-n)] = -1;
                }
                else{
                    m[y_bis - (i-n)][x - (j-n)] = input.get(x, y_bis);
                }
            }
        }

        return m;
    }

    public static void updatePixelNormalized(GrayU8 input, int m[][], int x, int y, int n){
        int sum = 0;
        int cpt = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                if(m[i][j] == -1){
                    cpt++;
                }
                else{
                    sum += m[i][j];
                }
            }
        }

        int nb_pixel = ((2*n) + 1) * ((2*n) + 1);
        nb_pixel -= cpt;
        sum /= nb_pixel;
        input.set(x, y, sum);
    }

    public static void meanFilterSimpleNormalized(Planar<GrayU8> image, Planar<GrayU8> imgout, int size) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int n = (size - 1)/2;

                for(int i = 0; i < input.height; i++) {
                    for(int j = 0; j < input.width; j++) {
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        int m[][] = new int[input.height][input.width];
                        m = updateMatNormalized(input,i,j,n);
                        updatePixelNormalized(output,m,j,i,n);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static int[][] updateMatExtended(GrayU8 input, int i, int j, int n){
        //y ligne x colonne
        int m[][] = new int[(2*n)+1][(2*n)+1];
        int y_tmp = 0;
        int x_tmp = 0;

        for(int y = i-n; y <= i+n; y++){
            for(int x = j-n; x <= j+n; x++){
                if(x < 0 || y < 0 || y >= input.height || x >= input.width){
                    if(x < 0){
                        x_tmp = 0;
                    }
                    if(y < 0){
                        y_tmp = 0;
                    }
                    if(y >= input.height){
                        y_tmp = input.height - 1;
                    }
                    if(x >= input.width){
                        x_tmp = input.width - 1;
                    }
                    m[y - (i-n)][x - (j-n)] = input.get(x_tmp, y_tmp);
                }
                else{
                    m[y - (i-n)][x - (j-n)] = input.get(x, y);
                }
            }
        }
        return m;
    }

    public static void updatePixelExtended(GrayU8 input, int m[][], int x, int y, int n){
        int sum = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                sum += m[i][j];
            }
        }

        int nb_pixel = ((2*n) + 1) * ((2*n) + 1);
        sum /= nb_pixel;
        //System.out.println(sum);
        input.set(x, y, sum);
    }

    public static void meanFilterSimpleExtended(Planar<GrayU8> image, Planar<GrayU8> imgout, int size) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int n = (size - 1)/2;

                for(int i = 0; i < input.height; i++) {
                    for(int j = 0; j < input.width; j++) {
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        int m[][] = updateMatExtended(input,i,j,n);
                        updatePixelExtended(output,m,j,i,n);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static int[][] updateMatReflect(GrayU8 input, int i, int j, int n){
        //y ligne x colonne
        int m[][] = new int[(2*n)+1][(2*n)+1];
        int y_tmp = 0;
        int x_tmp = 0;

        for(int y = i-n; y <= i+n; y++){
            for(int x = j-n; x <= j+n; x++){
                if(x < 0 || y < 0 || y >= input.height || x >= input.width){
                    x_tmp = x;
                    y_tmp = y;

                    if(x < 0){
                        x_tmp = -x;
                    }
                    if(y < 0){
                        y_tmp = -y;
                    }
                    if(y >= input.height){
                        y_tmp = (2*(input.height-1)) - y;
                    }
                    if(x >= input.width){
                        x_tmp = (2*(input.width-1)) - x;
                    }
                    m[y - (i-n)][x - (j-n)] = input.get(x_tmp, y_tmp);
                }
                else{
                    m[y - (i-n)][x - (j-n)] = input.get(x, y);
                }
            }
            }
        return m;
    }

    public static void updatePixelReflect(GrayU8 input, int m[][], int x, int y, int n){
        int sum = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                sum += m[i][j];
            }
        }

        int nb_pixel = ((2*n) + 1) * ((2*n) + 1);
        sum /= nb_pixel;
        //System.out.println(sum);
        input.set(x, y, sum);
    }

    public static void meanFilterSimpleReflect(Planar<GrayU8> image, Planar<GrayU8> imgout, int size) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int n = (size - 1)/2;

                for(int i = 0; i < input.height; i++) {
                    for(int j = 0; j < input.width; j++) {
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        int m[][] = updateMatReflect(input,i,j,n);
                        updatePixelReflect(output,m,j,i,n);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static void meanFilterWithBorders(Planar<GrayU8> input, Planar<GrayU8> output, int size, BorderType borderType) throws Exception{
        switch(borderType){
        case SKIP:
            meanFilterSimple(input,output,size);
            break;
        case EXTENDED:
            meanFilterSimpleExtended(input,output,size);
            break;
        case NORMALIZED:
            meanFilterSimpleNormalized(input,output,size);
            break;
        case REFLECT:
            meanFilterSimpleReflect(input,output,size);
            break;
        case WRAP:
            //NOT ASKED
            break;
        case ZERO:
            //NOT ASKED
            break;
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static void updatePixelConvolution(GrayU8 input, int m[][], int kernel[][], int x, int y, int n){
        int sum = 0;
        int cpt = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                sum += m[i][j] * kernel[i][j];
                cpt += kernel[i][j];
            }
        }

        int nb_pixel = cpt;
        sum /= nb_pixel;
        input.set(x, y, sum);
    }

    public static void convolution(Planar<GrayU8> image, Planar<GrayU8> imgout, int[][] kernel) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int size = kernel.length;
                int n = (size - 1)/2;

                for(int i = n; i < input.height - n; i++) {
                    for(int j = n; j < input.width - n; j++) {
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        int m[][];
                        m = updateMat(input,i,j,n);
                        updatePixelConvolution(output,m,kernel,j,i,n);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void gaussianFilter(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int kernel[][] = {
            {1,2,3,2,1},
            {2,6,8,6,2},
            {3,8,10,8,3},
            {2,6,8,6,2},
            {1,2,3,2,1}
        };

        convolution(image,imgout,kernel);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static int calculPixelGradient(int m[][], int kernel[][], int n){

        int sum = 0;

        for(int i = 0; i < (2*n) + 1; i++){
            for(int j = 0; j < (2*n) + 1; j++){
                sum += m[i][j] * kernel[i][j];
            }
        }

        return sum;
    }

    public static int[][] updatePixelsGradient(GrayU8 input, int[][] kernel) {
        int size = kernel.length;
        int n = (size - 1)/2;
        int res[][] = new int[input.height][input.width];

        for(int i = n; i < input.height - n; i++) {
            for(int j = n; j < input.width - n; j++) {
                int m[][] = updateMat(input,i,j,n);
                res[i][j] = calculPixelGradient(m,kernel,n);
            }
        }
        return res;
    }

    public static void gradientImageSobelColor(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception{
        int h1[][] = {
        {-1,0,1},
        {-2,0,2},
        {-1,0,1}
        };

        int h2[][] = {
        {-1,-2,-1},
        {0,0,0},
        {1,2,1}
        };

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int[][] res_h1;
                res_h1 = updatePixelsGradient(input,h1);

                int[][] res_h2;
                res_h2 = updatePixelsGradient(input,h2);

                int[][] res = new int[input.height][input.width];

                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        res[i][j] = (int)Math.sqrt( (res_h1[i][j]*res_h1[i][j]) + (res_h2[i][j]*res_h2[i][j]) );
                        output.set(j, i, res[i][j]);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void gradientImageSobelGrayLevel(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        Planar<GrayU8> tmp = image.createSameShape();

        colorToGray(image,tmp);

        int h1[][] = {
        {-1,0,1},
        {-2,0,2},
        {-1,0,1}
        };

        int h2[][] = {
        {-1,-2,-1},
        {0,0,0},
        {1,2,1}
        };

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = tmp.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                int[][] res_h1;
                res_h1 = updatePixelsGradient(input,h1);

                int[][] res_h2;
                res_h2 = updatePixelsGradient(input,h2);

                int[][] res = new int[input.height][input.width];

                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        res[i][j] = (int)Math.sqrt( (res_h1[i][j]*res_h1[i][j]) + (res_h2[i][j]*res_h2[i][j]) );
                        output.set(j, i, res[i][j]);
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////TP5//////////////////////////////////////////////////

    public static void colorToGray(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int[][] tmp = new int[image.height][image.width];

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                GrayU8 input = image.getBand(k_bis);

                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 0){ //RED
                            tmp[i][j] += (int)(0.3*input.get(j,i));
                            continue;
                        }
                        if(k_bis == 1){ //GREEN
                            tmp[i][j] += (int)(0.59*input.get(j,i));
                            continue;
                        }
                        if(k_bis == 2){ //BLUE
                            tmp[i][j] += (int)(0.11*input.get(j,i));
                            continue;
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor2.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 3){
                            output.set(j,i,input.get(j,i));
                            continue;
                        }

                        output.set(j,i,tmp[i][j]);
                    }
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void colorToGrayUnique(Planar<GrayU8> image, GrayU8 imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int[][] r = new int[image.height][image.width];
        int[][] g = new int[image.height][image.width];
        int[][] b = new int[image.height][image.width];

        for (int k = 0; k < image.getNumBands(); ++k){
            if(k >= 3){
                break;
            }
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);

                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 0){ //RED
                            r[i][j] += (int)(0.3*input.get(j,i));
                            continue;
                        }
                        if(k_bis == 1){ //GREEN
                            g[i][j] += (int)(0.59*input.get(j,i));
                            continue;
                        }
                        if(k_bis == 2){ //BLUE
                            b[i][j] += (int)(0.11*input.get(j,i));
                            continue;
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor2.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    imgout.set(j,i_bis,(int)(r[i_bis][j]*0.3 + g[i_bis][j]*0.59 + b[i_bis][j]*0.11));
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    static void rgbToHsv(int r, int g, int b, float[] hsv){
        float min = Math.min(Math.min(r,g), b);
        float max = Math.max(Math.max(r,g), b);

        for(int i=0;i<hsv.length;i++){
            switch(i){
                case 0: 
                    if(max==min){
                        hsv[i]=0;
                        break;
                    }
                    if(max==r){
                        hsv[i]=(60*(g-b)/(max-min)+360)%360;
                        break;
                    }
                    if(max==g){
                        hsv[i]=60*(b-r)/(max-min)+120;
                        break;
                    }
                    if(max==b){
                        hsv[i]=60*(r-g)/(max-min)+240;
                        break;
                    }
                    break;
                case(1):
                    if(max==0){
                        hsv[i]=0;
                        break;
                    }
                    else{
                        hsv[i]=1-(min/max);
                        break;
                    }
                case (2):
                    hsv[i]=max;
                    break;
                default:
                    break;
            }
        }
    }

    static void hsvToRgb(float h, float s, float v, int[] rgb){
        int hi=(int)(h/60)%6;
        float f = (float)(h/60)-hi;
        float l = (float)(v*(1-s));
        float m=(float)(v*(1-f*s));
        float n=(float)(v*(1-(1-f)*s));

        switch (hi){
            case 0:
                rgb[0]=(int)v;
                rgb[1]=(int)n;
                rgb[2]=(int)l;
                break;
            case 1:

                rgb[0]=(int)m;
                rgb[1]=(int)v;
                rgb[2]=(int)l;
                break;
            case 2:
                rgb[0]=(int)l;
                rgb[1]=(int)v;
                rgb[2]=(int)n;
                break;
            case 3:
                rgb[0]=(int)l;
                rgb[1]=(int)m;
                rgb[2]=(int)v;
                break;
            case 4:
                rgb[0]=(int)n;
                rgb[1]=(int)l;
                rgb[2]=(int)v;
                break;
            case 5:
                rgb[0]=(int)v;
                rgb[1]=(int)l;
                rgb[2]=(int)m;
                break;
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static void coloredFilter(Planar<GrayU8> image, Planar<GrayU8> imgout, int color) throws Exception{
        int num=image.getNumBands();
        int r[][] = new int[image.height][image.width];
        int g[][] = new int[image.height][image.width];
        int b[][] = new int[image.height][image.width];
        int a[][] = new int[image.height][image.width];

        int rgba[] = new int[3];
        float hsv[] = new float[3];


        if(num == 4){
            rgba = new int[4];
        }

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                
                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 0){ //RED
                            r[i][j] = input.get(j,i);
                        }
                        if(k_bis == 1){ //GREEN
                            g[i][j] = input.get(j,i);
                        }
                        if(k_bis == 2){ //BLUE
                            b[i][j] = input.get(j,i);
                        }
                        if(k_bis == 3){ //ALPHA
                            a[i][j] = input.get(j,i);
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        for (int k = 0; k < image.getNumBands(); ++k){
            final GrayU8 input = image.getBand(k);
            final GrayU8 output = imgout.getBand(k);
            
            for(int i = 0; i < input.height; i++){
                for(int j = 0; j < input.width; j++){
                    rgbToHsv(r[i][j],g[i][j],b[i][j],hsv);
                    if(num == 4){
                        rgba[3] = a[i][j]; 
                    }
                    hsvToRgb(color,hsv[1],hsv[2],rgba);
                    output.set(j,i,rgba[k]);
                }
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////

    public static int[] histogram(GrayU8 input){
        int histo[] = new int[256];

        for (int y = 0; y < input.height; ++y) {
            for (int x = 0; x < input.width; ++x) {
                int gl = input.get(x, y);
                histo[gl]++;
            }
        }

        return histo;
    }

    public static int[] histogramCumulated(GrayU8 input){
        int histo[] = histogram(input);
        int histoCumul[] = new int[256];

        for(int k = 0; k < 256; k++){
            for(int i = 0; i <= k; i++){
                histoCumul[k] += histo[i];
            }
        }

        return histoCumul;
    }

    public static void egalisationHistogram(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        Planar<GrayU8> tmp = image.createSameShape();

        colorToGray(image,tmp);

        int C[] = histogramCumulated(tmp.getBand(0));
        int N = tmp.height * tmp.width;

        for (int i = 0; i < image.getNumBands(); ++i){
            final int i_bis = i;
            executor.submit(() -> {
                final GrayU8 output = imgout.getBand(i_bis);
                final GrayU8 tmpImg = image.getBand(i_bis);
                for (int y = 0; y < output.height; ++y) {
                    for (int x = 0; x < output.width; ++x) {
                        int gl = tmpImg.get(x, y);
                        if(i_bis == 3){
                            output.set(x, y, gl);
                            continue;
                        }
                        else{
                            output.set(x, y, (int)(C[gl]*255/N));
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void egalisationHistogramUnique(GrayU8 image, GrayU8 imgout) throws Exception{
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int C[] = histogramCumulated(image);
        int N = image.height * image.width;
      
        for (int y = 0; y < imgout.height; ++y) {
            final int y_bis = y;
            executor.submit(() -> {
                for (int x = 0; x < imgout.width; ++x) {
                    int gl = image.get(x, y_bis);
                    imgout.set(x, y_bis, (int)(C[gl]*255/N));
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

////////////////////////////////////////////////////////////////////////////////////////////

    public static void egalisationHistogramV(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int r[][] = new int[image.height][image.width];
        int g[][] = new int[image.height][image.width];
        int b[][] = new int[image.height][image.width];

        float hsv[] = new float[3];

        int histo[] = new int[256];
        int histoCumul[] = new int[256];

        int N = image.height * image.width;

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                
                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 0){ //RED
                            r[i][j] = input.get(j,i);
                        }
                        if(k_bis == 1){ //GREEN
                            g[i][j] = input.get(j,i);
                        }
                        if(k_bis == 2){ //BLUE
                            b[i][j] = input.get(j,i);
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor2.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    rgbToHsv(r[i_bis][j],g[i_bis][j],b[i_bis][j],hsv);
                    histo[(int)hsv[2]]++;
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor3 = Executors.newFixedThreadPool(numThreads);

        for(int k = 0; k < 256; k++){
            final int k_bis = k;
            executor3.submit(() -> {
                for(int i = 0; i <= k_bis; i++){
                    histoCumul[k_bis] += histo[i];
                }
            });
        }

        executor3.shutdown();
        executor3.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor4 = Executors.newFixedThreadPool(numThreads);

        for (int i = 0; i < image.getNumBands(); ++i){
            final int i_bis = i;
            executor4.submit(() -> {
                final GrayU8 input = image.getBand(i_bis);
                final GrayU8 output = imgout.getBand(i_bis);
                for (int y = 0; y < input.height; ++y) {
                    for (int x = 0; x < input.width; ++x) {
                        int gl = input.get(x, y);
                        if(i_bis == 3){
                            output.set(x, y, gl);
                            continue;
                        }
                        else{
                            output.set(x, y, (int)(histoCumul[gl]*255/N));
                        }
                    }
                }
            });
        }

        executor4.shutdown();
        executor4.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }
      
//////////////////////////////////////////////////////////////////////////////////////////

    public static void egalisationHistogramS(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int r[][] = new int[image.height][image.width];
        int g[][] = new int[image.height][image.width];
        int b[][] = new int[image.height][image.width];

        float hsv[] = new float[3];

        int histo[] = new int[101];
        int histoCumul[] = new int[101];

        int N = image.height * image.width;

        for (int k = 0; k < image.getNumBands(); ++k){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                
                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis == 0){ //RED
                            r[i][j] = input.get(j,i);
                        }
                        if(k_bis == 1){ //GREEN
                            g[i][j] = input.get(j,i);
                        }
                        if(k_bis == 2){ //BLUE
                            b[i][j] = input.get(j,i);
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor2.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    rgbToHsv(r[i_bis][j],g[i_bis][j],b[i_bis][j],hsv);
                    histo[(int)(hsv[1]*100)]++;
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor3 = Executors.newFixedThreadPool(numThreads);

        for(int k = 0; k < 101; k++){
            final int k_bis = k;
            executor3.submit(() -> {
                for(int i = 0; i <= k_bis; i++){
                    histoCumul[k_bis] += histo[i];
                }
            });
        }

        executor3.shutdown();
        executor3.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor4 = Executors.newFixedThreadPool(numThreads);

        for (int i = 0; i < image.getNumBands(); ++i){
            final int i_bis = i;
            executor4.submit(() -> {
                final GrayU8 input = image.getBand(i_bis);
                final GrayU8 output = imgout.getBand(i_bis);
                for (int y = 0; y < input.height; ++y) {
                    for (int x = 0; x < input.width; ++x) {
                        int gl = input.get(x, y);
                        if(i_bis == 3){
                            output.set(x, y, gl);
                            continue;
                        }
                        else{
                            output.set(x, y, (int)(histoCumul[gl*100/256]*255/N));
                        }
                    }
                }
            });
        }

        executor4.shutdown();
        executor4.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static void negativeImage(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for (int k = 0; k < image.getNumBands(); ++k){
            final GrayU8 input = image.getBand(k);
            final GrayU8 output = imgout.getBand(k);
            final int k_bis = k;
            
            executor.submit(() -> {
            
                for(int i = 0; i < input.height; i++){
                    for(int j = 0; j < input.width; j++){
                        if(k_bis != 3){
                            output.set(j,i,255-input.get(j,i));
                        }
                        else{
                            output.set(j,i,input.get(j,i));
                        }                    
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

////////////////////////////////////////////////////////////////////////////////////////////

    public static int nodeSize(int k, int i, int j) throws Exception {
        String path = "./src/main/java/pdl/backend/AlgorithmImage/images/";
    
        if(k >= 0 && k <= 2500){
            path += "80x80/";
        }
        if(k > 2500 && k <= 10000){
            path += "40x40/";
        }
        if(k > 10000 && k <= 40000){
            path += "20x20/";
        }
        if(k > 40000 && k <= 197136){
            path += "9x9/";
        }
        if(k > 197136 && k <= 1776889){
            path += "3x3/";
        }
        if(k > 1776889){
            throw new Exception("The program encountered the following problem :\n\tk must be < 1776889\n");
        }
    
        //A for Alone, N for North, E for East, S for South, W for West.
    
        BufferedImage input = UtilImageIO.loadImage(path+"4.png");
        Planar<GrayU8> A = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        return A.width;
    }
  
    public static void createRandomGraph(int k, int m, int n, Planar<GrayU8> image) throws Exception{
        //Skin loading in terms of k :
        //0 to 2500 : 80x80 px
        //2501 to 10000 : 40x40 px
        //10001 to 40000 : 20x20 px
        //40001 to 197136 : 9x9 px
        //197137 to 1776889 : 3x3 px
        String path = "./src/main/java/pdl/backend/AlgorithmImage/images/";
    
        if(k >= 0 && k <= 2500){
            path += "80x80/";
        }
        if(k > 2500 && k <= 10000){
            path += "40x40/";
        }
        if(k > 10000 && k <= 40000){
            path += "20x20/";
        }
        if(k > 40000 && k <= 197136){
            path += "9x9/";
        }
        if(k > 197136 && k <= 1776889){
            path += "3x3/";
        }
        if(k > 1776889){
            throw new Exception("The program encountered the following problem :\n\tk must be < 1776889\n");
        }
    
        //A for Alone, N for North, E for East, S for South, W for West.
    
        BufferedImage input = UtilImageIO.loadImage(path+"0.png");
        Planar<GrayU8> NW = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"1.png");
        Planar<GrayU8> N = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"2.png");
        Planar<GrayU8> NE = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"3.png");
        Planar<GrayU8> W = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"4.png");
        Planar<GrayU8> A = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"5.png");
        Planar<GrayU8> E = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"6.png");
        Planar<GrayU8> SW = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"7.png");
        Planar<GrayU8> S = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        input = UtilImageIO.loadImage(path+"8.png");
        Planar<GrayU8> SE = ConvertBufferedImage.convertFromPlanar(input, null, true, GrayU8.class);
    
        if(k > m*n || m*A.height > 4000 || n*A.width > 4000){ //4000px = max size
            throw new Exception("The program encountered one of the following problems :\n\t- k must be <= m*n\n\t- m*(A.height) must be <= 4000\n\t- n*(A.width) must be <= 4000\n");
        }
    
        int[][] hashtable = new int[m][n];
        int cpt = 0;
    
        Random random = new Random();

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    
        for(int y = 0; y < m; y++){
            for(int x = 0; x < n; x++){
                if( (random.nextInt()%2) == 0 && cpt < k){
                    hashtable[y][x] = 16;
                    cpt++;
                }
                else{
                    hashtable[y][x] = 0;
                }
            }
        }
    
        while(cpt < k){
            for(int y = 0; y < m; y++){
                for(int x = 0; x < n; x++){
                    if(hashtable[y][x] == 0){
                        if( (random.nextInt()%2) == 0 ){
                            hashtable[y][x] = 16;
                            cpt++;
                        }
                    }
                }
            }
        }
    
        //Here we know that the hashtable is completely filled randomly.
        
        for(int y = 0; y < m; y++){
            final int y_bis = y;
            executor.submit(() -> {
                for(int x = 0; x < n; x++){
            
                    if(hashtable[y_bis][x] != 0){
            
                        if( (int)(x-1) >= 0 && (int)(y_bis-1) >= 0){ //NW
                            if(hashtable[(int)(y_bis-1)][(int)(x-1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 1;
                            }
                        }
                        if((int)(y_bis-1) >= 0){ //N
                            if(hashtable[(int)(y_bis-1)][x] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 2;
                            }
                        }
                        if((int)(x+1) < n && (int)(y_bis-1) >= 0){ //NE
                            if(hashtable[(int)(y_bis-1)][(int)(x+1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 4;
                            }
                        }
                        if((int)(x-1) >= 0){ //W
                            if(hashtable[y_bis][(int)(x-1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 8;
                            }
                        }
                        if((int)(x+1) < n ){ //E
                            if(hashtable[y_bis][(int)(x+1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 32;
                            }
                        }
                        if((int)(x-1) >= 0 && (int)(y_bis+1) < m){ //SW
                            if(hashtable[(int)(y_bis+1)][(int)(x-1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 64;
                            }
                        }
                        if((int)(y_bis+1) < m){ //S
                            if(hashtable[(int)(y_bis+1)][x] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 128;
                            }
                        }
                        if((int)(x+1) < n && (int)(y_bis+1) < m){ //SE
                            if(hashtable[(int)(y_bis+1)][(int)(x+1)] >= 16){
                                hashtable[y_bis][x] = hashtable[y_bis][x] | 256;
                            }
                        }
            
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    
        //Here final hashtable values are calculated.
    
        int square_height = A.height;
        int square_width = A.width;
    
        int r[][][] = new int[9][A.height][A.width];
        int g[][][] = new int[9][A.height][A.width];
        int b[][][] = new int[9][A.height][A.width];
        int a[][][] = new int[9][A.height][A.width];
    
        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        for(int w = 0; w < 9; ++w){
            final int w_bis = w;
            executor2.submit(() -> {
                for(int i = 0; i < A.height; ++i){
                    for(int j = 0; j < A.width; ++j){
                        for(int z = 0; z < A.getNumBands(); ++z){
                            GrayU8 in;

                            switch(w_bis){
                                case 0:
                                    in = NW.getBand(z);
                                    break;
                                case 1:
                                    in = N.getBand(z);
                                    break;
                                case 2:
                                    in = NE.getBand(z);
                                    break;
                                case 3:
                                    in = W.getBand(z);
                                    break;
                                case 4:
                                    in = A.getBand(z);
                                    break;
                                case 5:
                                    in = E.getBand(z);
                                    break;
                                case 6:
                                    in = SW.getBand(z);
                                    break;
                                case 7:
                                    in = S.getBand(z);
                                    break;
                                case 8:
                                    in = SE.getBand(z);
                                    break;
                                default:
                                    in = A.getBand(z);
                                    break;
                            }
                
                            switch(z){
                                case 0:
                                    r[w_bis][i][j] = in.get(j,i);
                                    break;
                                case 1:
                                    g[w_bis][i][j] = in.get(j,i);
                                    break;
                                case 2:
                                    b[w_bis][i][j] = in.get(j,i);
                                    break;
                                case 3:
                                    a[w_bis][i][j] = in.get(j,i);
                                    break;
                            }
                        }
                    }
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor3 = Executors.newFixedThreadPool(numThreads);
    
        for(int z = 0; z < image.getNumBands(); ++z){
            final int z_bis = z;
            executor3.submit(() -> {
                final GrayU8 output = image.getBand(z_bis);
            
                for(int cpt_y = 0; cpt_y < m; cpt_y++){
                    for(int cpt_x = 0; cpt_x < n; cpt_x++){
                
                        //NW
                        if( (hashtable[cpt_y][cpt_x] & 1) == 1){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[0][i][j] != 0 ||
                                    g[0][i][j] != 0 ||
                                    b[0][i][j] != 0 ||
                                    a[0][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[0][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[0][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[0][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[0][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //N
                        if( (hashtable[cpt_y][cpt_x] & 2) == 2){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[1][i][j] != 0 ||
                                    g[1][i][j] != 0 ||
                                    b[1][i][j] != 0 ||
                                    a[1][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[1][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[1][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[1][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[1][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //NE
                        if( (hashtable[cpt_y][cpt_x] & 4) == 4){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[2][i][j] != 0 ||
                                    g[2][i][j] != 0 ||
                                    b[2][i][j] != 0 ||
                                    a[2][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[2][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[2][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[2][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[2][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //W
                        if( (hashtable[cpt_y][cpt_x] & 8) == 8){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[3][i][j] != 0 ||
                                    g[3][i][j] != 0 ||
                                    b[3][i][j] != 0 ||
                                    a[3][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[3][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[3][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[3][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[3][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //A
                        if( (hashtable[cpt_y][cpt_x] & 16) == 16){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[4][i][j] != 0 ||
                                    g[4][i][j] != 0 ||
                                    b[4][i][j] != 0 ||
                                    a[4][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[4][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[4][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[4][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[4][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //E
                        if( (hashtable[cpt_y][cpt_x] & 32) == 32){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[5][i][j] != 0 ||
                                    g[5][i][j] != 0 ||
                                    b[5][i][j] != 0 ||
                                    a[5][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[5][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[5][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[5][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[5][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //SW
                        if( (hashtable[cpt_y][cpt_x] & 64) == 64){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[6][i][j] != 0 ||
                                    g[6][i][j] != 0 ||
                                    b[6][i][j] != 0 ||
                                    a[6][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[6][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[6][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[6][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[6][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //S
                        if( (hashtable[cpt_y][cpt_x] & 128) == 128){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[7][i][j] != 0 ||
                                    g[7][i][j] != 0 ||
                                    b[7][i][j] != 0 ||
                                    a[7][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[7][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[7][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[7][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[7][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        //SE
                        if( (hashtable[cpt_y][cpt_x] & 256) == 256){
                    
                            for(int i = 0; i < square_height; i++) {
                                for(int j = 0; j < square_width; j++){
                        
                                    if(r[8][i][j] != 0 ||
                                    g[8][i][j] != 0 ||
                                    b[8][i][j] != 0 ||
                                    a[8][i][j] != 0){
                                        switch(z_bis){
                                            case 0:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,r[8][i][j]);
                                                break;
                                            case 1:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,g[8][i][j]);
                                                break;
                                            case 2:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,b[8][i][j]);
                                                break;
                                            case 3:
                                                output.set(cpt_x*square_width+j,cpt_y*square_height+i,a[8][i][j]);
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }

        executor3.shutdown();
        executor3.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static int[][] createRandomGraphArray(int k, int m, int n) throws Exception{
        int[][] hashtable = new int[m][n];
        int cpt = 0;

        Random random = new Random();

        for(int y = 0; y < m; y++){
            for(int x = 0; x < n; x++){
                if( (random.nextInt()%(n*m)) == 0 && cpt < k){
                    hashtable[y][x] = 1;
                    cpt++;
                }
                else{
                    hashtable[y][x] = 0;
                }
            }
        }

        while(cpt < k){
            for(int y = 0; y < m; y++){
                for(int x = 0; x < n; x++){
                    if(hashtable[y][x] == 0){
                        if( (random.nextInt()%(n*m)) == 0 ){
                            hashtable[y][x] = 1;
                            cpt++;
                        }
                    }
                }
            }
        }

        return hashtable;
    }

    public static void zoneFilter(int k, int m, int n, Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception{
        if(k > m*n || m > image.height || n > image.width){
            throw new Exception("The program encountered one of the following problems :\n\tk must be <= m*n\n\tm must be <= image.height\n\tn must be <= image.width");
        }

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int[][] graph = createRandomGraphArray(k,m,n);
        float[][] graphS = new float[m][n];
        int [][] graphV = new int[m][n];

        int[][] r = new int[image.height][image.width];
        int[][] g = new int[image.height][image.width];
        int[][] b = new int[image.height][image.width];
        int[][] a = new int[image.height][image.width];

        int rgba[] = new int[3];
        float hsv[] = new float[3];

        int num = image.getNumBands();

        if(num == 4){
            rgba = new int[4];
        }
        //Obtain all RGB(A) values
        for(int l = 0; l < image.getNumBands(); l++){
            final int l_bis = l;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(l_bis);

                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        switch(l_bis){
                            case 0:
                                r[i][j] = input.get(j,i);
                                break;
                            case 1:
                                g[i][j] = input.get(j,i);
                                break;
                            case 2:
                                b[i][j] = input.get(j,i);
                                break;
                            case 3:
                                a[i][j] = input.get(j,i);
                                break;
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        //Fill the graph with values where a vertex is
        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor2.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    if(graph[i_bis*m/image.height][j*n/image.width] == 0){
                        continue;
                    }

                    rgbToHsv(r[i_bis][j],g[i_bis][j],b[i_bis][j],hsv);
                    graph[i_bis*m/image.height][j*n/image.width] += (int)hsv[0];
                    graphS[i_bis*m/image.height][j*n/image.width] += hsv[1];
                    graphV[i_bis*m/image.height][j*n/image.width] += (int)hsv[2];
                }
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor3 = Executors.newFixedThreadPool(numThreads);

        //Fill the graph
        for(int i = 0; i < m; i++){
            final int i_bis = i;
            executor3.submit(() -> {
                for(int j = 0; j < n; j++){
                    if(graph[i_bis][j] == 0){
                        int tmp_i = i_bis;
                        int tmp_j = j;

                        while(graph[tmp_i][tmp_j] == 0){
                            tmp_j++;
                            if(tmp_j%n != tmp_j){
                                tmp_i++;
                                tmp_i %= m;
                            }
                            tmp_j %= n;
                        }
                        graph[i_bis][j] = graph[tmp_i][tmp_j];
                        graphS[i_bis][j] = graphS[tmp_i][tmp_j];
                        graphV[i_bis][j] = graphV[tmp_i][tmp_j];
                    }
                }
            });
        }

        executor3.shutdown();
        executor3.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        //Fill the output image
        for (int l = 0; l < image.getNumBands(); l++){
            final GrayU8 output = imgout.getBand(l);
            
            for(int i = 0; i < image.height; i++){
                for(int j = 0; j < image.width; j++){
                    rgbToHsv(r[i][j],g[i][j],b[i][j],hsv);
                    if(num == 4){
                        rgba[3] = a[i][j]; 
                    }
                    hsv[0] = graph[i*m/image.height][j*n/image.width]/(image.height/m * image.width/n);
                    hsv[1] = graphS[i*m/image.height][j*n/image.width]/(image.height/m * image.width/n);
                    hsv[2] = graphV[i*m/image.height][j*n/image.width]/(image.height/m * image.width/n);
                    hsvToRgb(hsv[0],hsv[1],hsv[2],rgba);
                    output.set(j,i,rgba[l]);
                }
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////

    public static Boolean eraseElt(int tab[], int size, int elt) {

        for(int i = 0; i < size; i++){
            if(tab[i] == elt){
                for(int j = i; j < size-1; j++){
                    final int j_bis = j;
                    tab[j_bis] = tab[j_bis+1];
                }

                return true;
            }
        }
    
        return false;
    }
  
    public static void ShuffleGrid(Planar<GrayU8> image, Planar<GrayU8> imgout, int m, int n) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();

        Random random = new Random();
    
        int tab[] = new int[n*m];
        int size = n*m;
    
        int res[] = new int[n*m];
        int cpt = 0;
    
        int tmp, x, y;
    
        for(int i = 0; i < n*m; i++){
            final int i_bis = i;
            tab[i_bis] = i_bis;
        }
    
        while(size!=0){
            tmp = random.nextInt() % size;
            if(tmp < 0){
                tmp = -tmp;
            }
        
            //System.out.println("tmp= " + tmp + " (n*m)-size= " + (n*m-size) + "\n");
        
            res[(n*m)-size] = tab[tmp];
            if(eraseElt(tab,size,tab[tmp])){
                size--;
            }
        }
    
        int square_height = image.height/m;
        int square_width = image.width/n;
    
        for(int cpt_y = 0; cpt_y < m; cpt_y++){
            for(int cpt_x = 0; cpt_x < n; cpt_x++){
                tmp = res[cpt];
                x = tmp%n;
                y = (tmp-x)/n;
        
                ExecutorService executor = Executors.newFixedThreadPool(numThreads);

                for(int k = 0; k < image.getNumBands(); k++){
                    final int k_bis = k;
                    final int x_bis = x;
                    final int y_bis = y;
                    final int cpt_x_bis = cpt_x;
                    final int cpt_y_bis = cpt_y;
                    executor.submit(() -> {
                        final GrayU8 input = image.getBand(k_bis);
                        final GrayU8 output = imgout.getBand(k_bis);
                        for(int i = 0; i < square_height; i++) {
                            for(int j = 0; j < square_width; j++){
                                output.set(x_bis*square_width+j,y_bis*square_height+i, input.get(cpt_x_bis*square_width+j,cpt_y_bis*square_height+i));
                            }
                        }
                    });
                }

                executor.shutdown();
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        
                cpt++;
            }
        }
    }
  
//////////////////////////////////////////////////////////////////////////////////////////
  
    public static void ExtendImg(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    
        for(int k = 0; k < image.getNumBands(); k++){
            final int k_bis = k;
            executor.submit(() -> {
                int value0, value1,  value2, value3;

                final GrayU8 input = image.getBand(k_bis); 
                final GrayU8 output = imgout.getBand(k_bis);
            
                for (int x = 0; x < imgout.width; x++) {
                    for (int y = 0; y < imgout.height; y++) {
                
                        if(x%2 == 0 && y%2 == 0){
                            value0 = input.get(x/2, y/2);
                            output.set(x,y,value0);
                        }
                        else if( (x+y) % 2 == 0){
                            value0 = input.get( (x-1)/2, (y-1)/2);
                            value1 = input.get( (x+1)/2, (y-1)/2);
                            value2 = input.get( (x-1)/2, (y+1)/2);
                            value3 = input.get( (x+1)/2, (y+1)/2);
                            output.set(x,y, (value0+value1+value2+value3)/4 );
                        }
                        else{
                            
                            if(x > y){
                                value0 = input.get( (x-1)/2, y/2);
                                value1 = input.get( (x+1)/2, y/2);
                                output.set(x,y, (value0+value1)/2 );
                            }
                            else{
                                value0 = input.get(x/2, (y-1)/2);
                                value1 = input.get(x/2, (y+1)/2);
                                output.set(x,y, (value0+value1)/2 );
                            }
                        }
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    public static void horizontalMiror(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for(int k = 0; k < image.getNumBands(); k++){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        output.set(image.width-1 - j,i,input.get(j,i));
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void verticalMiror(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for(int k = 0; k < image.getNumBands(); k++){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        output.set(j,image.height-1 - i,input.get(j,i));
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

    public static void horizontalVerticalMiror(Planar<GrayU8> image, Planar<GrayU8> imgout) throws Exception {
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        for(int k = 0; k < image.getNumBands(); k++){
            final int k_bis = k;
            executor.submit(() -> {
                final GrayU8 input = image.getBand(k_bis);
                final GrayU8 output = imgout.getBand(k_bis);

                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        output.set(image.width-1 - j,image.height-1 - i,input.get(j,i));
                    }
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    public static boolean belongsTo(ArrayList<Integer> list, int n){
        for(int i=0; i<list.size(); i++){
            if(list.get(i) == n){
                return true;
            }
        }
        return false;
    }

    public static void LPE_updateNeigh(int[][] labels, int[][] stabilisation, int i, int j, int order){
        if(labels[i][j] == 0 || stabilisation[i][j] == 1){ //0 = WATERSHED
            return;
        }

        stabilisation[i][j] = 1; //FINAL STATE

        //Initialisation of neighbours list
        Point[] neighs = new Point[8];
        for(int k = 0; k < neighs.length; k++){
            final int k_bis = k;
            neighs[k_bis] = new Point();
        }    

        Point.updateNeighs(neighs,i,j);
        for(int k = 0; k < neighs.length; k++){
            //Case : the neighbour does not exist or has already a label
            if( neighs[k].getI() < 0 || neighs[k].getJ() < 0 || neighs[k].getI() >= labels.length || neighs[k].getJ() >= labels[0].length || 
            stabilisation[neighs[k].getI()][neighs[k].getJ()] != 0){ //0 = WATERSHED
                continue;
            }
            labels[neighs[k].getI()][neighs[k].getJ()] = labels[i][j];
            stabilisation[neighs[k].getI()][neighs[k].getJ()] = order;
        }
    }

    public static void LPE_homemade(Planar<GrayU8> image, Planar<GrayU8> imgout, int toler) throws Exception {
        if(toler < 0){
            throw new Exception("The program encountered the following problem :\n\ttolerance must be >= 0\n");
        }

        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        Planar<GrayU8> tmp = image.createSameShape();

        GrayU8 tmp1 = tmp.getBand(0).createSameShape();
        GrayU8 tmp2 = tmp.getBand(0).createSameShape();
        colorToGrayUnique(image,tmp1);
        egalisationHistogramUnique(tmp1,tmp2);

        int tolerance = toler; //20 is the "best" result giver integer

        int[][] labels = new int[image.height][image.width];
        int[][] stabilisation = new int[image.height][image.width]; //0 for "TO DO" & 1 for "DONE"

        //Initialisation of labels
        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    labels[i_bis][j] = -1; //INIT
                    stabilisation[i_bis][j] = 0; //TO DO
                }
            });
        }

        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor2 = Executors.newFixedThreadPool(numThreads);

        //Initialisation of the first label on the first pixel
        int label = 1;
        labels[0][0] = label; 

        //Initialisation of neighbours list
        Point[] neighs = new Point[8];
        for(int i = 0; i < neighs.length; i++){
            final int i_bis = i;
            executor2.submit(() -> {
                neighs[i_bis] = new Point();
            });
        }

        executor2.shutdown();
        executor2.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ArrayList<Integer> allLabels = new ArrayList<Integer>();
        ArrayList<ArrayList<Integer>> labelsValues = new ArrayList<ArrayList<Integer>>(); //RGB values
        ArrayList<Integer> labelsValuesGray = new ArrayList<Integer>(); //Gray values
        ArrayList<Integer> rgb0 = new ArrayList<Integer>();
        ArrayList<Integer> rgb1 = new ArrayList<Integer>();

        allLabels.add(0); //WATERSHED
        rgb0.add(255); rgb0.add(0); rgb0.add(0); //BLACK COLOR
        labelsValues.add(rgb0);
        labelsValuesGray.add(0);

        rgb1.clear();
        allLabels.add(1); //FIRST LABEL
        rgb1.add(image.getBand(0).get(0,0)); rgb1.add(image.getBand(1).get(0,0)); rgb1.add(image.getBand(2).get(0,0)); //FIRST PIXEL VALUES
        labelsValues.add(rgb1);
        labelsValuesGray.add(tmp2.get(0,0));

        //Labelise the entire image
        for(int i = 0; i < image.height; i++){
            for(int j = 0; j < image.width; j++){
                Point.updateNeighs(neighs,i,j);

                for(int k = 0; k < neighs.length; k++){
                    //Case : the neighbour does not exist or has already a label
                    if( neighs[k].getI() < 0 || neighs[k].getJ() < 0 || neighs[k].getI() >= image.height || neighs[k].getJ() >= image.width || 
                    labels[neighs[k].getI()][neighs[k].getJ()] == 0 ){ //0 = WATERSHED
                        continue;
                    }
                    //Same zone, same label
                    if( labelsValuesGray.get(labels[i][j]) - tolerance <= tmp2.get(neighs[k].getJ(),neighs[k].getI()) &&
                        tmp2.get(neighs[k].getJ(),neighs[k].getI()) <= labelsValuesGray.get(labels[i][j]) + tolerance &&
                        labels[i][j] != 0 //WATERSHED
                    ){
                        labels[neighs[k].getI()][neighs[k].getJ()] = labels[i][j];
                        continue;
                    }
                    //If the center pixel is a WATERSHED then give a label for each neighbour
                    if(labels[i][j] == 0){
                        label++;
                        labels[neighs[k].getI()][neighs[k].getJ()] = label;

                        //Determination of labels and its values
                        ArrayList<Integer> rgbn = new ArrayList<Integer>();
                        allLabels.add(label); //label
                        rgbn.clear();
                        rgbn.add(image.getBand(0).get(neighs[k].getJ(),neighs[k].getI())); 
                        rgbn.add(image.getBand(1).get(neighs[k].getJ(),neighs[k].getI())); 
                        rgbn.add(image.getBand(2).get(neighs[k].getJ(),neighs[k].getI()));
                        labelsValues.add(rgbn); //DETERMINE A UNIQUE COLOR FOR THE LABEL
                        labelsValuesGray.add(tmp2.get(neighs[k].getJ(),neighs[k].getI()));
                        continue;
                    }
                    //Else it is a watershed
                    labels[neighs[k].getI()][neighs[k].getJ()] = 0;
                }
            }
        }
        //End of labelising

        //Stabilisation of the output
        int cpt = 1;

        //Stabilisation of the watersheds
        while(cpt != 0){
            cpt = 0;

            //Update and complete the watershed
            for(int i=0; i<image.height; i++){
                for(int j=0; j<image.width; j++){

                    if(labels[i][j] == 0){ //Already a part of the watershed
                        continue;
                    }

                    int nonExistingNeigh = 0;
                    int nbNeighWatershed = 0;
                    int nbCardNeighWatershed = 0;
                    Point.updateNeighs(neighs,i,j);

                    for(int k = 0; k < neighs.length; k++){

                        if( neighs[k].getI() >= 0 && neighs[k].getJ() >= 0 && 
                        neighs[k].getI() < image.height && neighs[k].getJ() < image.width 
                        ){ 

                            if(labels[neighs[k].getI()][neighs[k].getJ()] == 0){ //0 = WATERSHED
                                nbNeighWatershed++;

                                if(Point.isCard(neighs[k],i,j)){
                                    nbCardNeighWatershed++;
                                }
                            }
                        }
                        else{ //The neighbour does not exist
                            nonExistingNeigh++;
                        }

                        if(nbNeighWatershed >= 6 || nbCardNeighWatershed >= 3 || (nbNeighWatershed >= 2 && nonExistingNeigh >= 3) ){ //The last : useful to handle minor errors on non-watershed pixels (which would be watersheded) 
                            labels[i][j] = 0; //Become a part of the watershed
                            cpt++;
                            break;
                        }

                    }
                }
            }
        }

        ExecutorService executor3 = Executors.newFixedThreadPool(numThreads);

        //Initialisation of cardinal neighbours list
        Point[] cardNeighs = new Point[4];
        for(int i = 0; i < cardNeighs.length; i++){
            final int i_bis = i;
            executor3.submit(() -> {
                cardNeighs[i_bis] = new Point();
            });
        }

        executor3.shutdown();
        executor3.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor4 = Executors.newFixedThreadPool(numThreads);

        //Stabilisation of the bassins
        for(int i = 0; i < image.height; i++){
            final int i_bis = i;
            executor4.submit(() -> {
                for(int j = 0; j < image.width; j++){
                    if(labels[i_bis][j] == 0){ //WATERSHED
                        stabilisation[i_bis][j] = 1; //WATERSHED = STABILISED
                    }
                }
            });
        }

        executor4.shutdown();
        executor4.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        ExecutorService executor5 = Executors.newFixedThreadPool(numThreads);

        cpt = 1;
        int cpt2 = 1;
        int order = 2;

        while(cpt != 0){ //cpt = nb of unstabilised + non-treated cases left
            cpt = 0;

            while(cpt2 != 0){ //cpt2 = nb of stabilised but next cases to be used left (considered as non-treated)
                cpt2 = 0;

                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        if(stabilisation[i][j] == order){ //order (case which will change)
                            LPE_updateNeigh(labels,stabilisation,i,j,order+1); //Determination of next orders cases
                            cpt2++;
                            cpt++;
                        }
                    }
                }

                order++;
            }

            for(int i = 0; i < image.height; i++){
                if(cpt != 0){
                    break;
                }

                for(int j = 0; j < image.width; j++){
                    if(stabilisation[i][j] == 0){ //NOT STABILISED YET
                        LPE_updateNeigh(labels,stabilisation,i,j,order); 
                        cpt++;
                        cpt2++; //To enter in the second while
                        break;
                    }
                }
            }
        }

        //Colorisation of the output image
        for(int k = 0; k < image.getNumBands(); k++){
            final int k_bis = k;
            executor5.submit(() -> {
                final GrayU8 output = imgout.getBand(k_bis);
                for(int i = 0; i < image.height; i++){
                    for(int j = 0; j < image.width; j++){
                        if(k_bis == 3){ //A in RGBA is always equal to 100% (1)
                            output.set(j,i,255);
                        }
                        else{
                            output.set( j,i,labelsValues.get(labels[i][j]).get(k_bis) );
                        }
                    }
                }
            });
        }

        executor5.shutdown();
        executor5.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    public static void PopArt(Planar<GrayU8> input, Planar<GrayU8> output,int backgroundColor) {
     
        ThresholdImageOps.threshold(input.getBand(0), output.getBand(0),(int) (GThresholdImageOps.computeOtsu(input.getBand(0), 0, 255)), true);

        for (int y = 0; y < output.getHeight(); y++) {
            for (int x = 0; x < output.getWidth(); x++) {
                int rgb = output.getBand(0).get(x, y);
                if (rgb == 0) { // Background color
                    rgb = backgroundColor;
                }
              
                float[] hsv =new float[3];
                ColorHsv.rgbToHsv((rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF, hsv);
                hsv[0] *= 2; // Increase the hue
                hsv[1] *= 2; // Increase the saturation
                int rgbs = ColorHsv.hsvToRgb(hsv[0], hsv[1], hsv[2]);
                for (int b = 0; b < output.getNumBands(); b++) {
                    output.getBand(b).set(x, y, rgbs & 0xFF);
                    rgbs >>= 8;
                }
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static void main( String[] args ) throws Exception {   
        //Load image
		if (args.length < 2) {
			System.out.println("missing input or output image filename");
			System.exit(-1);
		}
		final String inputPath = args[0];

        BufferedImage image = UtilImageIO.loadImage(inputPath);
        Planar<GrayU8> input = ConvertBufferedImage.convertFromPlanar(image, null, true, GrayU8.class);
        Planar<GrayU8> output = input.createSameShape();

        //Call here a function
        egalisationHistogramV(input, output);

        final String outputPath = args[1];
		UtilImageIO.saveImage(output, outputPath);
		System.out.println("Image saved in: " + outputPath);
    }
}
