#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#define X_SIZE 1280
#define Y_SIZE 720 
#define SLIDERS 6
#define CX_SIZE 100*SLIDERS
#define CY_SIZE 400
#define STOP 20
#define SBOT CY_SIZE-150
#define S1LOC CY_SIZE-100
#define S2LOC CY_SIZE-50

/* here are our X variables */
Display *display;
int screen;
float bf,rf,gf,PI;
Window win,cwin;
GC gc,cgc;
XImage *x_image;
XImage *cx_image;
unsigned char *x_buffer;
unsigned char *cx_buffer;
unsigned char *buff;
unsigned char *buff2;
unsigned long black,white;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();
void *lfo();
void *dbuf();
void *control();
void slider(int);
void knob(int,int,int);
void toggle(int,int,int);
float osc (float, float, float , float, int ,int);
float ra,ga,ba;
float byf,ryf,gyf;
int rw,gw,bw;
int rph,gph,bph;
int rw1,gw1,bw1;
int rmax,gmax,bmax,ready;
int tick,rate;
int state[2];
int mix,del;


/* Jpegs */
int read_JPEG_file (char *, unsigned char *, int *);
int jayit(unsigned char *,int, int, char *);

void usage ()
{
	printf("usage: synth\n");
	exit (1);
}

int main(int argc,char *argv[])
{
	PI=2*M_PI;
        struct timespec tim, tim2;
        tim.tv_sec = 0;
       	tim.tv_nsec = 10000000;

        buff=(unsigned char *)malloc(sizeof(unsigned char)*6*X_SIZE*Y_SIZE);
	pthread_t lfo_id,dis_id,cnt_id;
	pthread_create(&lfo_id, NULL, lfo, NULL);
	pthread_create(&dis_id, NULL, dbuf, NULL);
	pthread_create(&cnt_id, NULL, control, NULL);
	init_x();
	char dum[10];
	int XS,TS;
	int x,y,r;
	XS=X_SIZE*3;
	TS=XS*Y_SIZE;
	bf=1; rf=1; gf=1; ryf=1; gyf=1; byf=1; ba=127; ga=127; ra=127; rate=2; r=0; state[0]=0; state[1]=0; rph=0; gph=0; bph=0; mix=50; del=399;
	while (1>0)
	{
		while (state[0]>0 && state[1]>0)
		{
       			nanosleep(&tim , &tim2);
		}
		if (state[0]==0 && state[1]==0)
		{
			r=0;
		}
		if (state[0]==1 && state[1]==0)
		{
			r=1;
			state[r]=2;
			//printf ("two %d \n",r);
		}
		if (state[0]==0 && state[1]==1)
		{
			r=0;
			state[r]=2;
			//printf ("two %d \n",r);
		}
		for (y=0;y<Y_SIZE;y++)
		{
			int Y;
			float bbf,rrf,ggf,rp,gp,bp;
			bp=osc(0,X_SIZE,bph,y,Y_SIZE,0);
			rp=osc(0,X_SIZE,rph,y,Y_SIZE,0);
			gp=osc(0,X_SIZE,gph,y,Y_SIZE,0);
			bbf=osc(0,bmax,byf,y,Y_SIZE,bw);
			rrf=osc(0,rmax,ryf,y,Y_SIZE,rw);
			ggf=osc(0,gmax,gyf,y,Y_SIZE,gw);
			Y=y*XS+(r*TS);
			for (x=0;x<X_SIZE;x++)
			{
				int p;
				p=Y+(3*x);
				buff[p]=osc(0,ra,rrf,x+rp,X_SIZE,rw);
				buff[p+1]=osc(0,ga,ggf,x+gp,X_SIZE,bw);
				buff[p+2]=osc(0,ba,bbf,x+bp,X_SIZE,gw);
			}
		}
		if ( state[0]==2 && state[1]==0){ state[0]=1;}
		else if ( state[1]==2 && state[0]==0){ state[1]=1;}
		else { state[r]=1; }
		
	}
}

float osc (float min, float max, float cycles, float p, int span,int form)
{
	float ret,diff;
	diff=max-min;
	if (form==0)
	{
		ret=min+(diff*(1+sin(PI*p*cycles/span))/2);
	}
	else if (form==1)
	{
		ret=min+((diff*((int)(p*cycles)%span))/span);
	}
	else if (form==2)
	{
		if (((float)(((int)(p*cycles))%span)/span)<0.5){ ret=min;}else{ret=max;}
	}else if (form==3)
	{
		ret=min+rand()%(int)diff;
	}

	return ret;
}

void *dbuf ()
{
	struct timespec tim, tim2;
	int m,in,r,wri,rea,n;
	unsigned char *verb;
	m=X_SIZE*Y_SIZE*4;
	n=X_SIZE*Y_SIZE*3;
	tim.tv_sec = 0;
        tim.tv_nsec = 10000000L;
	wri=del;
	rea=0;
        verb=(unsigned char *)malloc(sizeof(unsigned char)*400*n);
       	while (1>0)
	{
		in=tick;
       		// wait 10ms for ready
		while ( state[0]==0 && state[1]==0)
		{
       			nanosleep(&tim , &tim2);
		}
		if (state[0]==1){ r=0;}
		if (state[1]==1){ r=1;}
		int p,q,rd,wr,pix;
		q=r*n;
		rd=rea*n;
		wr=wri*n;
		pix=100-mix;
		for (p=0;p<m;p+=4)
		{
			char nr,ng,nb;
			nr=((pix*verb[rd++])+(mix*buff[q++]))/100;
			x_buffer[p+2]=nr;
			ng=((pix*verb[rd++])+(mix*buff[q++]))/100;
			x_buffer[p+1]=ng;
			nb=((pix*verb[rd++])+(mix*buff[q++]))/100;
			x_buffer[p]=nb;
			verb[wr++]=nb;
			verb[wr++]=nr;
			verb[wr++]=ng;
		}
		rea++;if( rea>del){rea=0;}	
		wri=rea+del;if(wri>del){wri=wri-del-1;}
		//if(wri>0){rea=wri-1;}else{rea=99;}
		//memcopy(buff
		//printf ("rea %d wri %d del %d\n",rea,wri,del);
		state[r]=0;
		while (tick-in<rate)
		{
       			nanosleep(&tim , &tim2);
		}
		XPutImage(display, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	}
}

void *lfo ()
{
       struct timespec tim, tim2;
       tim.tv_sec = 0;
       tim.tv_nsec = 10000000L;
       nanosleep(&tim , &tim2);

	while (1>0)
	{
                nanosleep(&tim , &tim2);
		tick++;
		ba=osc(10,255,1,tick,1403,0);
		ra=osc(10,255,1,tick,1070,0);
		ga=osc(18,255,1,tick,1100,0);

		rph=osc(0,5,1,tick,2000,1);
		gph=osc(0,9,1,tick,3400,1);
		bph=osc(0,13,1,tick,2900,1);

		gw=osc(0,3,1,tick,4200,1);
		bw=osc(0,3,1,tick,3100,1);
		rw=osc(0,3,1,tick,5000,1);

		gyf=osc(0,3,1,tick,30400,1);
		byf=osc(0,3,1,tick,30410,1);
		ryf=osc(0,3,1,tick,30420,1);

		bmax=osc(0,50,1,tick,9060,0);
		gmax=osc(0,50,1,tick,8000,0);
		rmax=osc(0,50,1,tick,7200,0);

		//mix=osc(20,30,1,tick,3000,0);
		//del=osc(399,1,1,tick,6110,1);
		//del=18;
		//mix=40;

		//scanf("%d",&bf);
		//scanf("%d",&rf);
		//scanf("%d",&gf);

	}
}

void *control()
{
	XEvent event;           /* the XEvent declaration !!! */
	struct timespec tim, tim2;
       	tim.tv_sec = 2;
       	tim.tv_nsec = 10000000L;
	int x_point,y_point,kpoint,ppoint1,ppoint2,sl;
	
        nanosleep(&tim , &tim2);
	printf("Starting  \n");
	kpoint=STOP;
	ppoint1=STOP;
	ppoint2=STOP;
	for (sl=0;sl<SLIDERS;sl++){ slider(sl);knob(kpoint,sl,1);}
	while (1>0)
	{
             	XNextEvent(display, &event);
                if (event.type==ButtonPress ) {
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			printf("press x %d y %d \n",x_point,y_point);
		}
                if (event.type==ButtonRelease ) {
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			printf("release x %d y %d \n",x_point,y_point);
			
		}
		if (event.type=MotionNotify){
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			//printf("move  x %d y %d \n",x_point,y_point);
			if (y_point>10 && y_point<CY_SIZE-10)
			{
				if (x_point<CX_SIZE/2){
					knob(ppoint1,0,0);
					knob(y_point,0,1);
					ppoint1=y_point;
					mix=100*(y_point-10)/(CY_SIZE-20);
					printf("mix %d\n",mix);
				}else{
					knob(ppoint2,0,0);
					knob(y_point,0,1);
					ppoint2=y_point;
					del=300*(y_point-10)/(CY_SIZE-20);
					printf("delay %d\n",del);
				}
			}
		}
	}
}

void slider(int num)
{
	int loop,xp;
	xp=(CX_SIZE*(num+1))/(SLIDERS+1);
	for (loop=STOP;loop<SBOT;loop++)
	{
		XDrawPoint(display, cwin, cgc, xp, loop);
	}
}

void knob(int val,int num,int m)
{
	int xp;
	if (m==1){XSetForeground(display,cgc,white);}else{XSetForeground(display,cgc,black);}
	xp=(CX_SIZE*(num+1))/(SLIDERS+1);
	int h,w,hc,wc;
	h=4;w=14;
	for (hc=val-h;hc<=val+h;hc++)
	{
		for (wc=3;wc<=w;wc++)
		{
			XDrawPoint(display, cwin, cgc, xp-wc, hc);
			XDrawPoint(display, cwin, cgc, xp+wc, hc);
		}
	}
}

void toggle(int val,int num,int m)
{
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
	XInitThreads();
        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        cx_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*CX_SIZE*CY_SIZE);
        display=XOpenDisplay((char *)0);
        screen=DefaultScreen(display);
        black=BlackPixel(display,screen),
        white=WhitePixel(display,screen);
        win=XCreateSimpleWindow(display,DefaultRootWindow(display),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        cwin=XCreateSimpleWindow(display,DefaultRootWindow(display),0,0,
                CX_SIZE, CY_SIZE, 5, white,black);
        XSetStandardProperties(display,win,"image","images",None,NULL,0,NULL);
        XSetStandardProperties(display,cwin,"control","control",None,NULL,0,NULL);
        XSelectInput(display, cwin, ExposureMask|ButtonPressMask|KeyPressMask|ButtonReleaseMask|ButtonMotionMask);
        gc=XCreateGC(display, win, 0,0);
        cgc=XCreateGC(display, cwin, 0,0);
        XSetBackground(display,gc,black); XSetForeground(display,gc,white);
        XSetBackground(display,cgc,black); XSetForeground(display,cgc,white);
        XClearWindow(display, win);
        XClearWindow(display, cwin);
        XMapRaised(display, win);
        XMapRaised(display, cwin);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(display, 0);
        x_image=XCreateImage(display, visual, DefaultDepth(display,DefaultScreen(display)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
        //cx_image=XCreateImage(display, visual, DefaultDepth(display,DefaultScreen(display)), ZPixmap, 0, cx_buffer, CX_SIZE, CY_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(display, gc);
        XDestroyWindow(display,win);
        XCloseDisplay(display);
        exit(1);
};

void redraw() {
        XClearWindow(display, win);
};


