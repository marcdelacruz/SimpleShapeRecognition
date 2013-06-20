
/**
 * Author: Marc DelaCruz
 * Project 2 Computer Vision CS 6384 Summer 2011
 */
import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javax.swing.*;
import javax.media.jai.*;
import com.sun.media.jai.widget.*;
import java.math.*;
import java.util.*;
import java.lang.*;

public class proj2 extends JFrame {
    public static boolean SHOW_DEBUGS = false;
    public proj2(String img1Name, String img2Name) {
        
	// create the JFrame with a title
	super(img1Name + " " + img2Name);
	addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		    System.exit(0);
		}
	    });
	getContentPane().setLayout(new FlowLayout());
        PlanarImage img1 = JAI.create("fileload", img1Name);
	getContentPane().add(new DisplayJAI(img1));
        PlanarImage img2 = JAI.create("fileload", img2Name);
	getContentPane().add(new DisplayJAI(img2));

	// Input cannot be color image
	if(img1.getSampleModel().getNumBands() != 1
	   || img2.getSampleModel().getNumBands() != 1) {
	    System.out.println("Input should be gray-level images");
	    System.exit(0);
	}
	// convert the image to an array of integers
	int[][] IMG1 = imageBandTOarray(img1,0);
	int[][] IMG2 = imageBandTOarray(img2,0);

        pack();
        setVisible(true);

	double[] info1 = ComputeInformation(IMG1);
	double[] info2 = ComputeInformation(IMG2);
	AnalyzeInformation(info1, info2);
    }

    double[] ComputeInformation(int[][] img) {
	// The information computed by this method is the image moments.
	// You may want to change it to compute additional information
	// if you believe it is necessary.

	double m00=0.0, m01=0.0, m10=0.0, m20=0.0, m11=0.0, m02=0.0;

	int width = img[0].length;
	int height = img.length;
    
        
	for(int i = 0 ; i < height ; i++)
	    for(int j = 0 ; j < width ; j++) {
        int p = img[i][j] > 0 ? 1 : 0;
		m00 += p;
		m10 += j*p; 
		m01 += i*p;
		m20 += j*j*p;
		m11 += j*i*p;
		m02 += i*i*p;
        
	    }
        
        
        double M00 = m00;
        double M10 = 0;
        double M01 = 0;
        double M11 = m11 - ((m10*m01)/m00);
        double M20 = m20 - ((m10*m10)/m00);
        double M02 = m02 - ((m01*m01)/m00);
        
        double scaleFactor = (Math.sqrt(m00));

        double I1 = M20 + M02;
        double I2 = (Math.pow(M20 - M02 , 2.0)) + 4*Math.pow(M11,2.0);
        
        double M20g = (I1 + Math.sqrt(I2))/(2.0*Math.pow(M00,2.0));
        
        double M02g = (I1 - Math.sqrt(I2))/(2.0*Math.pow(M00,2.0));
        
        final double TRIANGLE_M20_M02G = Math.round(0.096*1000);
        final double CIRCLE_M20_M02G = Math.round(0.0795*1000);
        double isCircle = (CIRCLE_M20_M02G == Math.round(M20g*1000) &&
                           CIRCLE_M20_M02G == Math.round(M20g*1000)) ? 1 : 0;
        
        double isTriangle = (TRIANGLE_M20_M02G == Math.floor(M20g*1000) &&
                             TRIANGLE_M20_M02G == Math.floor(M20g*1000)) ? 1 : 0;
        double height_ = 0;
        double width_ = 0;
        double slant = 0; //0 = horizontal, 1 = slant, 2 = vertical
        if (SHOW_DEBUGS)
        {
            System.out.println("~~~~~~~~M20=" + M20 + " M02=" + M02 + " M11=" + M11);
            System.out.println("~~~~~~~~M20g=" + M20g + " M02g=" + M02g + " Math.floor(M20g*1000)=" + Math.floor(M20g*1000));
        }
        double _2Alpha = (M20 - M02) != 0.0 ? -(2.0*M11)/(M20 - M02) : 0.0;
        double angle = Math.atan(_2Alpha)/2.0;

        if (angle == Double.NaN || ((M20 - M02) == 0.0 && isCircle == 0) )
        {
            //multiply by the sin 0f _2Alpha
            angle = (22.0/28.0);
            if (SHOW_DEBUGS)
            {
                System.out.println(">>>~~~~~~~~~~~~~~~~~~~angle=" + angle*57.2957795 +
                                   " degrees _2alpha=" + _2Alpha + " (Math.sin(_2Alpha)=" + (adjustedSinRad(_2Alpha)));
            }
        }
        double rot_angle = angle;

        final double TO_DEGREES = 57.2957795;

        double angleDegrees = angle*TO_DEGREES;
        System.out.println("width=" + width + " height=" + height);
        System.out.println("m00=" + m00);
        System.out.println("m10=" + m10 + "  m01=" + m01);
        System.out.println("m20=" + m20 + "  m11=" + m11 + "  m02=" + m02);

        if (SHOW_DEBUGS)
        {
            System.out.println("~~~~~~~~~~~~~~~~~~~angle=" + angleDegrees +
                               " degrees _2alpha=" + _2Alpha);
        }
        int counter = 0;
        boolean hasFoundAngle = false;
        int[] smoothingList = new int[9];
        do
        {
            if (SHOW_DEBUGS)
            {
                System.out.println("counter=" + counter);
            }
            angle = angle + (isTriangle == 0 ? ((22.0/14.0)*(counter > 0 ? 1 : 0)) : ((22.0/14.0)*(counter > 0 ? 1 : 0)));
            if (SHOW_DEBUGS)
            {
                System.out.println("~~~~~~~~Angle=" + angle*TO_DEGREES + "(added=" + ((22.0/14.0)*counter)*TO_DEGREES + ")");
            }
            int[][] rotatedImage = new int[height][width];

            for(int i = 0 ; i < height ; i++)
            {
                for(int j = 0 ; j < width ; j++) {
                    int x0 = (int)Math.round(m10/m00);
                    int y0 = (int)Math.round(m01/m00);
/*
                    int X = (int)Math.round( (((j-x0)*adjustedCosRad(angle)) + 
                                            ((i-y0)*adjustedSinRad(angle))  
                                              +x0
                                              ) );
                    int Y = (int)Math.round( ((-(j-x0)*adjustedSinRad(angle)) +
                                             ((i-y0)*adjustedCosRad(angle)) 
                                              +y0
                                              ) );
                    
                    if (X >= 0 && Y >= 0 && X < width && Y < height)
                    {
                        //rotatedImage[Y][X] = 1;
                        rotatedImage[i][j] = img[Y][X] > 0 ? 1 : 0;
                    }
                    */
                    
                    
                    //transformation
                     double x = (double)(((j-x0)*Math.cos(angle)) + 
                                          ((i-y0)*Math.sin(angle))  
                                             +x0
                                             ) ;
                     double y = (double)((-(j-x0)*Math.sin(angle)) +
                                         ((i-y0)*Math.cos(angle)) 
                                         +y0
                                         );

                    int x_flr = (int)Math.floor(x);
                    int y_flr = (int)Math.floor(y);
                    int x_ceil = (int)Math.ceil(x_flr);
                    int y_ceil = (int)Math.ceil(y_flr);
                    double a = (double)(x - x_flr);
                    double b = (double)(y - y_flr);
                    if (x_ceil < width && y_ceil < height && 
                        x_flr < width && y_flr < height &&
                        x_ceil >= 0 && y_ceil >= 0 && 
                        x_flr >= 0 && y_flr >= 0)
                    {
                    rotatedImage[i][j] = (int)Math.round (
                                            ((1-a)*(1-b)*(img[y_flr][x_flr ] > 0 ? 1 : 0)  ) +
                                            ((1-a)*(b)*(img[ y_ceil][ x_flr ] > 0 ? 1 : 0)  ) +
                                            ((a)*(1-b)*(img[ y_flr ][ x_ceil] > 0 ? 1 : 0)  ) +
                                            ((a)*(b)*(img[y_ceil][ x_ceil ]> 0 ? 1 : 0)  ) );
                    }
                    
                }
            }
            /*
            for(int i = 0 ; i < height ; i++) //apply smoothing
            {
                for(int j = 0 ; j < width ; j++)
                {
                    Arrays.fill(smoothingList, 0);
                    int count = 0;
                    int loop = 0;
                    for (int k = -1 ; k < 2; k++)
                    {
                        for (int l= -1; l < 2; l++)
                        {
                            if (i+k >= 0 && i+k < height &&
                                j+l >= 0 && j+l < width)
                            {
                                //add the number to list
                                smoothingList[count++] = rotatedImage[i+k][j+l];
                            }
                            else
                            {
                                //add zero
                                smoothingList[count++] = 0;
                            }
                        }
                    }

                    Arrays.sort(smoothingList);
                    rotatedImage[i][j] = smoothingList[4];
                }
            }*/
            
        double m00_=0.0, m01_=0.0, m10_=0.0, m20_=0.0, m11_=0.0, m02_=0.0;
        for(int i = 0 ; i < height ; i++)
            for(int j = 0 ; j < width ; j++) {
                int p = rotatedImage[i][j];
                m00_ += p;
                m10_ += j*p; 
                m01_ += i*p;
                m20_ += j*j*p;
                m11_ += j*i*p;
                m02_ += i*i*p;
            }
            for(int i = 0 ; i < height ; i++)
                for(int j = 0 ; j < width ; j++) {
                    rotatedImage[i][j] = rotatedImage[i][j] > 0 ? 255 : 0;
                }
            PlanarImage p = arrayTOimage(rotatedImage);
            getContentPane().add(new DisplayJAI(p));
            
        double _m20 = m20_ - ((m10_*m10_)/m00_);
        double _m02 = m02_ - ((m01_*m01_)/m00_);
            
            width_ = Math.pow((Math.pow(_m20,3.0)*144.0)/_m02,
                              (1.0/8.0));
            height_ = (12.0*_m20)/Math.pow(width_,3.0);
            
            // compute the angle again , it should be zero now if rotated correctly
            double _2Alpha2 = (_m20 - _m02) != 0.0 ? -(2.0*(m11_ - ((m10_*m01_)/m00_)))/(_m20 - _m02) : 0.0;
            double angle_ = Math.atan(_2Alpha2)/2.0;
            if (SHOW_DEBUGS)
            {
                System.out.println(">>>>> TRanslated angle_=" + angle_*TO_DEGREES);
            }
            if (Math.round(angle_*TO_DEGREES) == 0 || 
                (Math.abs(Math.round(angle_*TO_DEGREES)) < 5 )  )
            {
                hasFoundAngle = true;
                break;
            }
            
            counter++;
        }
        while ( counter < 5 && !hasFoundAngle);
        
        if ((Math.abs(angleDegrees)%90) <= 30 && (Math.abs(angleDegrees)%90) >= 0 )
            slant = 0;
        else if ((Math.abs(angleDegrees)%90) <= 60 && (Math.abs(angleDegrees)%90) >= 30 )
        {
            slant = 1;

            if (SHOW_DEBUGS)
            {
                System.out.println("~~~~~~~~~~~~~~~~ SLANT angle=" + Math.abs(angle*TO_DEGREES));
            }
        }
        rot_angle = (-rot_angle)+((Math.PI/2.0)*counter); //add 90 degrees
        if (width_ < height_ && width_ != height_)
        {
            double temp = width_;
            width_ = height_;
            height_ = temp;
            if (slant == 0 && (isTriangle == 0))
                slant = 2;
            rot_angle = rot_angle+(Math.PI/2.0); //add 90 degrees
            if (SHOW_DEBUGS)
            {
                System.out.println("~~~~~~~~~~~~~~~~ Added 90 degrees");
            }

        }
        else
        {
            
        }
        angleDegrees = angle*TO_DEGREES;

        double triangleSide = 
                   Math.round(
                     Math.sqrt( (4.0*(Math.floor(height_-1))*(Math.floor(width_-1)))/Math.sqrt(3.0) )-0.25
                              ) ;
        double circleRadius = Math.round(Math.sqrt( (7.0*height_*width_)/22.0) );
        if (SHOW_DEBUGS)
        {
            System.out.println("~~~~~~~~~~~~~~~~~~~height=" + height_ + 
                               " width=" + width_ + 
                               " triangleSide=" + triangleSide +
                               " circleRadius=" + circleRadius + " counter=" + counter);
        }

	// the following information need not be printed
	System.out.println("width=" + width + " height=" + height);
	System.out.println("m00=" + m00);
	System.out.println("m10=" + m10 + "  m01=" + m01);
	System.out.println("m20=" + m20 + "  m11=" + m11 + "  m02=" + m02);

        
	double[] info = //{width, height, m00, m01, m10, m20, m11, m02};
               {width_, height_, isCircle, circleRadius, 
                  isTriangle, triangleSide, slant, rot_angle, 
                 Math.floor(M20g*1000), Math.floor(M02g*1000),
                 Math.round(m10/m00), Math.round(m01/m00)
                 };

        //width height circle radius triangle side, slant, orientation, M20G, M02G, centerX, centerY
        
        
	return(info);
    }

    void AnalyzeInformation(double[] info1, double[] info2) {

	// Your code should be here
        final double TO_DEGREES = 57.2957795;
    /**********************************************************/
        
    /**********************************************************/
        String roughOrientationStr[] = { "HORIZONTAL", "SLANT", "VERTICAL" };
    // The following information should be printed:
	System.out.println
	    ("******* Information about the object in first image:");

        System.out.println("Object location is: (" + info1[10] + "," +
                                                     info1[11] +")");
	System.out.println
	    ("Rough object orientation: " + roughOrientationStr[ (int)info1[6] ]);
	System.out.println("Exact orientation: " + info1[7] + 
                       " radians (" + Math.round(info1[7]*TO_DEGREES) + " degrees)");
        
    System.out.println("The object " + (info1[2] == 1 ? "IS" : "ISN'T") + " a circle");
	System.out.println("If the object is a circle its radius is: " + info1[3]); 
	System.out.println("The object " + (info1[4] == 1 ? "IS" : "ISN'T") + " an equilateral triangle");
	System.out.println
	    ("If it is an equilateral triangle its side length is: " + info1[5]);

    System.out.println
     ("If it is a rectangle its width=" + info1[0] + " height=" + info1[1]);                                                        
	System.out.println
	    ("******* Information about the object in second image:");


    System.out.println("Object location is: (" + info2[10] + "," +
                      info2[11] +")");
	System.out.println
	    ("Rough object orientation:  " + roughOrientationStr[ (int)info2[6] ]);
        System.out.println("Exact orientation: " + info2[7] + 
                           " radians (" + Math.round(info2[7]*TO_DEGREES) + " degrees)");
    System.out.println("The object " + (info2[2] == 1 ? "IS" : "ISN'T") + " a circle");
    System.out.println("If the object is a circle its radius is: " + info2[3]); 
    System.out.println("The object " + (info2[4] == 1 ? "IS" : "ISN'T") + " an equilateral triangle");
    System.out.println
    ("If it is an equilateral triangle its side length is: " + info2[5]);
    System.out.println
    ("If it is a rectangle its width=" + info2[0] + " height=" + info2[1]); 

    boolean hasSameShape = (info1[8] == info2[8] && info1[9] == info2[9]) ? true : false;
	System.out.println
	    ("The objects in first image and the object in the second image "
	     +  (hasSameShape ? "HAVE" : "DO NOT HAVE") + " the same shape");
    }

    private int[][] imageBandTOarray(PlanarImage image, int band) {
	Raster raster = image.getData();
        int width = raster.getWidth();
        int height = raster.getHeight();
	int bands = raster.getNumBands();

	if(band > bands || band < 0) {
	    System.out.println("Can't handle band=" +band+ " bands=" +bands);
	    System.exit(0);
	}
	int[][] array = new int[height][width];
	for(int i = 0 ; i < height ; i++)
	    raster.getSamples(0 , i, width, 1, band, array[i]);
	return(array);
    }
    
    private PlanarImage arrayTOimage(int[][] array) {
        // Creates a gray level image. 
        // Values in array are assumed to be in 0-255 range
        int height = array.length;
        int width = array[0].length;
        BufferedImage bimg = new BufferedImage 
	    (width, height, BufferedImage.TYPE_BYTE_GRAY);
        WritableRaster raster = bimg.getRaster();
        for(int i = 0 ; i < height ; i++)
            raster.setSamples(0, i, width, 1, 0, array[i]);
        return(PlanarImage.wrapRenderedImage(bimg));
    }
    double adjustedSinRad(double angleInRad)
    {
        return adjustedSin(angleInRad*(180/Math.PI));
    }
    double adjustedCosRad(double angleInRad)
    {
        return adjustedCos(angleInRad*(180/Math.PI));
    }
    double adjustedSin(double angleInDegrees)
    {
        double sint = 0.0;
        if (angleInDegrees ==  -270) { sint = 1; }
        else if(angleInDegrees ==  -180) { sint = 0; }
        else if(angleInDegrees ==  -90) {  sint = -1;  }
        else if(angleInDegrees ==  0) {   sint = 0;  }
        else if(angleInDegrees ==  90) {  sint = 1;  }
        else if(angleInDegrees ==  180) { sint = 0;  }
        else if(angleInDegrees ==  270) { sint = -1;  }
        else if(angleInDegrees ==  360) { sint = 0;  }
        else {
            double tt = angleInDegrees*Math.PI/180;
            sint = Math.sin(tt);
        }
        return sint;
    }
    
    double adjustedCos(double angleInDegrees)
    {
        double cost = 0.0;
        if (angleInDegrees ==  -270) { cost = 0; }
        else if(angleInDegrees ==  -180) { cost = -1; }
        else if(angleInDegrees ==  -90) { cost = 0; }
        else if(angleInDegrees ==  0) {  cost = 1; }
        else if(angleInDegrees ==  90) { cost = 0; }
        else if(angleInDegrees ==  180) { cost = -1; }
        else if(angleInDegrees ==  270) { cost = 0; }
        else if(angleInDegrees ==  360) { cost = 1; }
        else {
            double tt = angleInDegrees*Math.PI/180;
            cost = Math.cos(tt);
        }
        return cost;
    }
    
    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage: proj2 Image1 Image2");
            System.exit(1);
        }
	new proj2(args[0], args[1]);
    }
}
