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
#define HIST 30
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

struct slid {
        float value[4];
        char name[8];
        int wv[3];
        int pos[4];
        int mcol;
        float smax[3];
        float smin[3];
        float max[3];
        float min[3];
        float osc[3];
        int eam[4];
};

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();
void *lfo();
void *dbuf();
void *control();
void slider(struct slid *,int);
void knob(int,int,int);
void c_toggle(struct slid *,int,int);
void b_toggle(struct slid *,int,int,int,int,int);
float osc (float, float, float , float, int ,int);
float param[10][3];
float byf,ryf,gyf;
int rw,gw,bw;
float rph,gph,bph;
float rmax,gmax,bmax;
int tick,rate;
int state[2];
int mix,del;
int *mask;
XColor    color[100];
struct slid *slids;

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

        slids=(struct slid *)malloc(sizeof(struct slid)*(SLIDERS+10));
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
	bf=1; rf=1; gf=1; ryf=1; gyf=1; byf=1; rate=2; r=0; state[0]=0; state[1]=0; rph=0; gph=0; bph=0; mix=50; del=HIST-1;
	while (1>0)
	{
		float amplitude[3];
		float phase[3];
		float wavform[3];
		float wavfrm[3];
		float freq[3];
		float freqc[3];
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
		float bbf,rrf,ggf,rp,gp,bp;
		// copy in the values so they dont change during the frame.
		int cl;
		for (cl=0;cl<3;cl++)
		{
			amplitude[cl]=param[2][cl];
			freq[cl]=param[3][cl];
			phase[cl]=param[4][cl];
			wavform[cl]=param[5][cl];
			wavfrm[cl]=param[6][cl];
			freqc[cl]=param[9][cl];
		}

		for (y=0;y<Y_SIZE;y++)
		{
			float maxf[3];
			for (cl=0;cl<3;cl++){ maxf[cl]=osc(0,freq[cl],freqc[cl],y,Y_SIZE,wavfrm[cl]); }
			int Y;
			Y=y*XS+(r*TS);
			for (x=0;x<X_SIZE;x++)
			{
				int p;
				p=Y+(3*x);
				for (cl=0;cl<3;cl++){ buff[p+cl]=osc(0,amplitude[cl],maxf[cl],x+phase[cl],X_SIZE,wavform[cl]); }
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
	// sin
	if (form==0)
	{
		ret=min+(diff*(1+sin(PI*p*cycles/span))/2);
	}
	else if (form==1)
	{
		ret=min+((diff*((int)(p*cycles)%span))/span);
	}
	// square
	else if (form==2)
	{
		if (((float)(((int)(p*cycles))%span)/span)<0.5){ ret=min;}else{ret=max;}
	// rand
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
        verb=(unsigned char *)malloc(sizeof(unsigned char)*(HIST+1)*n);
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
			verb[wr++]=nr;
			verb[wr++]=ng;
			verb[wr++]=nb;
		}
		rea++;if( rea>del){rea=0;}	
		wri=rea+del;if(wri>del){wri=wri-del-1;}
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
		/*ba=osc(10,255,1,tick,4403,0);
		ra=osc(10,255,1,tick,5070,0);
		ga=osc(18,255,1,tick,7100,0);*/
		// oscillate sliders 2-SLIDERS 1 is the mixer.
		int a;
		for (a=2;a<=SLIDERS;a++)
		{
			int clr;
			for (clr=0;clr<3;clr++)
			{
				if (slids[a].eam[clr]==1){ param[a][clr]=osc(slids[a].smin[clr],slids[a].smax[clr],1,tick,slids[a].osc[clr],0);}
			}
		}

		rph=osc(0,1,1,tick,2000,0);
		gph=osc(0,1,1,tick,2000,0);
		bph=osc(0,1,1,tick,2000,0);

		gw=osc(0,3,1,tick,9000,1);
		bw=osc(0,3,1,tick,9000,1);
		rw=osc(0,3,1,tick,9000,1);

		param[9][0]=osc(0,5,1,tick,3400,0);
		param[9][1]=osc(0,5,1,tick,3400,0);
		param[9][2]=osc(0,5,1,tick,3400,0);

		//bmax=osc(0,6.1,1,tick,9060,0);
	//	gmax=osc(0,7.2,1,tick,8000,0);
	//	rmax=osc(0,8.2,1,tick,7200,0);

		//mix=osc(20,30,1,tick,3000,0);
		//del=osc(HIST,1,1,tick,6110,1);
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
	int *bval,clr;
	struct timespec tim, tim2;
       	tim.tv_sec = 1;
       	tim.tv_nsec = 10000000L;
	int x_point,y_point,kpoint,ppoint1,ppoint2,sl,n;
	char eam[10];

	strcpy(eam,"extoscman");

        bval=(int *)malloc(sizeof(int)*(SLIDERS+10));
        mask=(int *)malloc(sizeof(int)*CX_SIZE*CY_SIZE);

	strcpy(slids[1].name,"mix  ");
	strcpy(slids[2].name,"amp  ");
	strcpy(slids[3].name,"freq ");
	strcpy(slids[4].name,"phase");
	strcpy(slids[5].name,"wvfm1");
	strcpy(slids[6].name,"wvfm2");

	for (clr=0;clr<3;clr++) { 
		slids[2].min[clr]=0; slids[2].max[clr]=255;  
		slids[3].min[clr]=0; slids[3].max[clr]=300; 
		slids[4].min[clr]=0; slids[4].max[clr]=X_SIZE;  
		slids[5].min[clr]=0; slids[5].max[clr]=3;  
		slids[6].min[clr]=0; slids[6].max[clr]=3; 
		slids[2].smin[clr]=slids[2].min[clr]; slids[2].smax[clr]=slids[2].max[clr];
		slids[3].smin[clr]=slids[3].min[clr]; slids[3].smax[clr]=slids[3].max[clr];
		slids[4].smin[clr]=slids[4].min[clr]; slids[4].smax[clr]=slids[4].max[clr];
		slids[5].smin[clr]=slids[5].min[clr]; slids[5].smax[clr]=slids[5].max[clr];
		slids[6].smin[clr]=slids[6].min[clr]; slids[6].smax[clr]=slids[6].max[clr];
	}

	for (n=0;n<CX_SIZE*CY_SIZE;n++){mask[n]=0;}
	
        nanosleep(&tim , &tim2);
	printf("Starting  \n");
	// BAsic layout
	for (sl=1;sl<=SLIDERS;sl++){ slider(slids,sl);knob(STOP,sl,1);
		c_toggle(slids+sl,sl,0);
		b_toggle(slids+sl,sl,sl*CX_SIZE/(SLIDERS+3),S2LOC,3,0);
		slids[sl].eam[0]=2;
		slids[sl].eam[1]=2;
		slids[sl].eam[2]=2;
		slids[sl].eam[3]=2;
		XSetForeground(display,cgc,white);
		XClearArea(display, cwin, ((sl*CX_SIZE)/(SLIDERS+3))+18,S2LOC-10,20,20,0 );
		XDrawString(display,cwin,cgc,((sl*CX_SIZE)/(SLIDERS+3))+18,S2LOC,eam+(3*(slids[sl].eam[slids[sl].mcol])),3);
	}
	b_toggle(slids+10,210,CX_SIZE-90,40,3,0);
	b_toggle(slids+11,211,CX_SIZE-90,80,3,0);
	int pressed;pressed=0;
	while (1>0)
	{
		int rval,t;
             	XNextEvent(display, &event);
                if (event.type==ButtonPress ) {
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			rval=mask[x_point+(y_point*CX_SIZE)];
			if (rval>0 && rval<100){pressed=1;}
			printf("press x %d y %d val %d\n",x_point,y_point,rval);
			if (rval>100 && rval<200) { t=rval-100; c_toggle(slids+t,t,1); }
			if (rval>200 ) { t=rval-200; 
				b_toggle(slids+t,t,t*CX_SIZE/(SLIDERS+3),S2LOC,2,1);}
		}
                if (event.type==ButtonRelease ) {
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			printf("release x %d y %d \n",x_point,y_point);
			if (rval>100 && rval<200) { c_toggle(slids+t,rval-100,0); }
			if (rval>200 ) { t=rval-200;b_toggle(slids+t,t,t*CX_SIZE/(SLIDERS+3),S2LOC,2,0);}
			printf ("t %d mcol %d eam %d\n",t,slids[t].mcol,slids[t].eam[slids[t].mcol]);
			XSetForeground(display,cgc,white);
			XClearArea(display, cwin, ((t*CX_SIZE)/(SLIDERS+3))+18,S2LOC-10,20,20,0 );
			XDrawString(display,cwin,cgc,((t*CX_SIZE)/(SLIDERS+3))+18,S2LOC,eam+(3*(slids[t].eam[slids[t].mcol])),3);
			pressed=0;
		}
		if (event.type==MotionNotify && pressed==1){
                        x_point=event.xbutton.x; y_point=event.xbutton.y;
			//printf("move  x %d y %d \n",x_point,y_point);
			if (y_point>=STOP && y_point<=SBOT)
			{
				knob(y_point,rval,1);
				if (slids[rval].eam[slids[rval].mcol]==1){slids[rval].osc[slids[rval].mcol]=y_point*10;}else{
				if (rval==1 && slids[1].mcol==0){mix=(100*(y_point-STOP))/(SBOT-STOP);}
				else if (rval==1 && slids[1].mcol==1){del=(HIST*(y_point-STOP))/(SBOT-STOP);}
				else {
					if (slids[rval].mcol!=3){param[rval][slids[rval].mcol]=(slids[rval].max[0]*(y_point-STOP))/(SBOT-STOP);}else{
					param[rval][0]=(slids[rval].max[0]*(y_point-STOP))/(SBOT-STOP);
					param[rval][1]=(slids[rval].max[1]*(y_point-STOP))/(SBOT-STOP);
					param[rval][2]=(slids[rval].max[2]*(y_point-STOP))/(SBOT-STOP);
					}
				}
				}
				printf("mix %d\n",mix);
				printf("del %d\n",del);
			}
		}
	}
	free (mask);
}

void slider(struct slid *slids,int num)
{
	int loop,xp;
	xp=(CX_SIZE*num)/(SLIDERS+3);
	XDrawString(display,cwin,cgc,xp-6,10,slids[num].name,5);
	for (loop=STOP-5;loop<SBOT+5;loop++)
	{
		XDrawPoint(display, cwin, cgc, xp, loop);
	}
}

void knob(int yval,int num,int m)
{
	int xp,x,y,h,w;
	h=4;w=14;
	XSetForeground(display,cgc,black);
	xp=(CX_SIZE*num)/(SLIDERS+3);
	for (y=STOP-h;y<=SBOT+h;y++)
	{
		int mp; mp=CX_SIZE*y;
		for (x=3;x<=w;x++)
		{
			int x1,x2; x1=xp-x;x2=xp+x;
			XDrawPoint(display, cwin, cgc, x1, y);
			XDrawPoint(display, cwin, cgc, x2, y);
			mask[mp+x1]=0; mask[mp+x2]=0;
		}
	}
	XSetForeground(display,cgc,white);
	for (y=yval-h;y<=yval+h;y++)
	{
		int mp; mp=CX_SIZE*y;
		//3 is the 3 pixels either side of the two slider knob parts
		for (x=3;x<=w;x++)
		{
			int x1,x2; x1=xp-x;x2=xp+x;
			XDrawPoint(display, cwin, cgc, x1, y);
			XDrawPoint(display, cwin, cgc, x2, y);
			mask[mp+x1]=num; mask[mp+x2]=num;
		}
	}
}

void c_toggle(struct slid *col,int num,int m)
{
	int xp;
	xp=(CX_SIZE*num)/(SLIDERS+3);
	int h,w,x,y,c;
	h=8;w=8;
	c=col->mcol;
	if (m==0){
		XSetForeground(display,cgc,color[c].pixel);
		for (y=-h;y<=h;y++)
		{
			int mp,yy; yy=y+S1LOC;mp=CX_SIZE*yy;
			for (x=-w;x<=w;x++)
			{
				XDrawPoint(display, cwin, cgc, xp+x, yy);
				mask[mp+x+xp]=num+100; 
			}
		}
	}else{
		col->mcol++;if(col->mcol>3){col->mcol=0;}
		XSetForeground(display,cgc,black);
		for (y=1-h;y<h;y++)
		{
			int yy; yy=y+S1LOC;
			for (x=1-w;x<w;x++)
			{
				XDrawPoint(display, cwin, cgc, xp+x, yy);
			}
		}
	}
}

void b_toggle(struct slid *slids,int num,int xp,int yp, int max,int m )
{
        int h,w,x,y;
        h=8;w=8;
        if (m==0){
                slids->eam[slids->mcol]++;if(slids->eam[slids->mcol]>max){slids->eam[slids->mcol]=0;}
                XSetForeground(display,cgc,white);
                for (y=-h;y<=h;y++)
                {
                        int mp,yy; yy=y+yp;mp=CX_SIZE*yy;
                        for (x=-w;x<=w;x++)
                        {
                                XDrawPoint(display, cwin, cgc, xp+x, yy);
                                mask[mp+x+xp]=num+200;
                        }
                }
        }else{
                XSetForeground(display,cgc,black);
                for (y=1-h;y<h;y++)
                {
                        for (x=1-w;x<w;x++)
                        {
                                XDrawPoint(display, cwin, cgc, xp+x, yp+y);
                        }
                }
	}
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
	//
	//
	//
        Colormap cmap;
        cmap = DefaultColormap(display, screen);
        color[3].red = 65535; color[3].green = 65535; color[3].blue = 65535;
        color[0].red = 65535; color[0].green = 0; color[0].blue = 0;
        color[1].red = 0; color[1].green = 65535; color[1].blue = 0;
        color[2].red = 0; color[2].green = 0; color[2].blue = 65535;
        XAllocColor(display, cmap, &color[0]);
        XAllocColor(display, cmap, &color[1]);
        XAllocColor(display, cmap, &color[2]);
        XAllocColor(display, cmap, &color[3]);

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


