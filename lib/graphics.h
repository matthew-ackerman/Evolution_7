#include <GL/freeglut.h>
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <GL/glx.h>			/* GLXContext */
#include <GL/freeglut_std.h>

#define NIL 		(0)
#define X_SIZE 		1060
#define Y_SIZE 		1000 
#define MAX_COL 	65536			/* Max number of colours.	*/
#define SCALE		100.0
#define CELL		1

#define	COLOR_A 20
#define COLOR_B 2
#define COLOR_C 0.2

#define ORTH	0
#define PERS	1

void init(void);
void display();

GLint window;           /* The number of our GLUT window */
GLint Xsize=800;
GLint Ysize=800;
GLfloat dist=20.0;

GLfloat xangle=0.0,yangle=0.0,zangle=0.0;   /* axis angles */

int LastX=-1;
int LastY=-1;
int mode=0;

double X_length=1.0, Y_length=1.0, Z_length=1.0, x_trans=0.0, y_trans=0.0;
bool rotate_flag=false,drag_x=false, drag_y=false, drag_z=false, translate_flag=false, x_proj=false, y_proj=false, z_proj=true;

tensor <double> *T;
LIST *list;

/* Simple window transformation routine */

GLvoid Transform(GLfloat Width, GLfloat Height){
	if(mode==PERS){
		glViewport(0, 0, Width, Height);              /* Set the viewport */
		glMatrixMode(GL_PROJECTION);                  /* Select the projection matrix */
		glLoadIdentity();				/* Reset The Projection Matrix */
		gluPerspective(90.0,Width/Height,0.1,100.0);  /* Calculate The Aspect Ratio Of The Window */
		glMatrixMode(GL_MODELVIEW);                   /* Switch back to the modelview matrix */
	}
	else if (mode==ORTH){
		glViewport (0, 0, Width, Height);
		glMatrixMode (GL_PROJECTION);
		glLoadIdentity ();
		glOrtho(-15.0, 15.0, -15.0, 15.0, -15.0, 15.0);
		gluPerspective (90.0, Width/Height, 0.1, 100.0);
		glMatrixMode (GL_MODELVIEW);
		gluOrtho2D(0, Width, 0, Height);
	}
}

void keyboard(unsigned char key, int x, int y){
	if (key==27) { //27 is the ascii code for the ESC key
		exit (0); //end the program
	}
	else if (key=='x'){
		xangle=-90.0;
		yangle=-90.0;
		x_proj=true;
		y_proj=false;
		z_proj=false;
	}
	else if (key=='y'){
		xangle=0;
		yangle=90.0;
		x_proj=false;
		y_proj=true;
		z_proj=false;
	}
	else if (key=='z'){
		xangle=0.0;
		yangle=0.0;
		x_proj=false;
		y_proj=false;
		z_proj=true;
	}
	else if (key=='1'){
		mode=ORTH;
		Transform(Xsize,Ysize);
	}
	else if (key=='2'){
		mode=PERS;
		Transform(Xsize,Ysize);
	}
}


void MouseFunc(int button, int state, int x, int y) {
	if (state==0){
		cout << "Press button:" << button << endl;
		if(button==0){
			LastX=x;
			LastY=y;
			rotate_flag=true;
		}
		if(button==2){
			LastX=x;
			LastY=y;
			if((x-Xsize/2.0)/(y-Ysize/2.0)< 1&& (x-Xsize/2.0)/(y-Ysize/2.0) > 0){
				drag_x=true;
			}
			if((x-Xsize/2.0)/(y-Ysize/2.0)< 0 && (x-Xsize/2.0)/(y-Ysize/2.0) > -1){
				drag_y=true;
			}
			if((x-Xsize/2.0)/(y-Ysize/2.0)< -2 && (x-Xsize/2.0)/(y-Ysize/2.0) > -3){
				drag_x=true;
			}	
		}
		if(button==3){
			dist+=(.1*dist);
			cout << dist << endl;
			cout << x_trans << endl;
		}
		else if(button==4){
			dist-=(.1*dist);
			cout << dist << endl;
			cout << x_trans << endl;
		}
		else if(button==1){
			LastX=x;
			LastY=y;
			translate_flag=true;
		}
	}
	if (state==1){
		if(button==0){
			rotate_flag=false;
		}
		if(button==2){
			drag_x=false;
			drag_y=false;
			drag_z=false;
		}
		if(button==1){
			translate_flag=false;
		}
	}
}

void MotionFunc(int x, int y) {

	if(rotate_flag){
		xangle+=(double)(x-LastX)/10.0;
		yangle+=(double)(y-LastY)/10.0;
	}
	if(drag_x){
		X_length+=(double)(x-LastX)/100.0;
	}
	if(drag_y){
		Y_length+=(double)(x-LastX)/100.0;
	}
	if(drag_z){
		Z_length+=(double)(x-LastX)/100.0;
	}
	if(translate_flag){
		x_trans+=(double)(x-LastX)/10.0;
		y_trans-=(double)(y-LastY)/10.0;		
	}

//	zangle+=(double)
	LastX=x;
	LastY=y;
	
}

/* A general OpenGL initialization function.  Sets all of the initial parameters. */
GLvoid InitGL(GLfloat Width, GLfloat Height) {
	
	Xsize=Width;
	Ysize=Height;
	glClearColor(1.0, 1.0, 1.0, 1.0);	/* This Will Clear The Background Color To Black */
	glLineWidth(2.0);                       /* Add line width,   ditto */
	Transform( Width, Height );             /* Perform the transformation */
//	glEnable(GL_DEPTH_TEST);
        glutPostRedisplay();
	glutMotionFunc(MotionFunc);
	glutMouseFunc(MouseFunc);
	glutKeyboardFunc(keyboard);
}

/* The function called when our window is resized  */
GLvoid ReSizeGLScene(GLint Width, GLint Height) {
  if (Height==0)    Height=1;                   /* Sanity checks */
  if (Width==0)      Width=1;
  Transform( Width, Height );                   /* Perform the transformation */
}

void setDrawTensor(tensor <double> * Set){
	T=Set;
}

void (*UpdateTensor) (void);

void setUpdateTensor(void ( *function ) (void)){
	UpdateTensor=*function;
}

GLvoid DrawGLSceneTensor() {
	UpdateTensor();


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	/* Clear The Screen And The Depth Buffer */

	glPushMatrix();
	glLoadIdentity();

	glTranslatef(0.0,0.0,-dist);
	glRotatef(xangle,0.0,1.0,0.0);
	glRotatef(yangle,1.0,0.0,0.0);

	INDEX index(2);
	double valueA, valueB, valueC, meanValue;

	double width=(double)(T->get_width(0));
	double height=(double)(T->get_width(1));

	double max;

	if (width>height) max=width;
	else max=height;

	double rwidth=width/max;
	double rheight=height/max;

#define WIDTH (rwidth)
#define HEIGHT (rheight)

	glBegin(GL_LINES);

	glColor3f(0.0,0.0,0.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(15.0*cos((xangle/180.0)*M_PI),15*sin((yangle/180.0)*M_PI),15.0*sin((((xangle)/180.0))*M_PI));

	glColor3f(1.0,0.0,0.0);
	glRotatef(xangle,0.0,1.0,0.0);
	glRotatef(yangle,1.0,0.0,0.0);

	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,15*sin((yangle/180.0)*M_PI),0.0);

	glColor3f(0.0,0.0,1.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(0.0,0.0,15.0*sin((((xangle)/180.0))*M_PI));

	glColor3f(0.0,1.0,1.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(15.0*cos((xangle/180.0)*M_PI),0.0,0.0);
	glEnd();

	glBegin(GL_LINES);

	glColor3f(1.0,0.0, 0.0);

	glVertex3f(15.0*(WIDTH)-((CELL-0.5)*30.0/height),0,15.0*(HEIGHT)-((CELL)*30.0/width));
	glVertex3f(15.0*(WIDTH)-((CELL-0.5)*30.0/height),0,-15.0*(HEIGHT));

	glColor3f(0.0,1.0, 0.0);

	glVertex3f(15.0*(WIDTH)-((CELL)*30.0/height),0,-15.0*(HEIGHT));
	glVertex3f(-15.0*(WIDTH),0,-15.0*(HEIGHT));

	glColor3f(0.0,0.0, 1.0);

	glVertex3f(-15.0*(WIDTH),0,-15.0*(HEIGHT));
	glVertex3f(-15.0*(WIDTH),0,15.0*(HEIGHT)-((CELL)*30.0/width));

	glColor3f(0.5,0.5, 0.5);

	glVertex3f(-15.0*(WIDTH),0,15.0*(HEIGHT)-((CELL)*30.0/width));
	glVertex3f(15.0*(WIDTH)-((CELL)*30.0/height),0,15.0*(HEIGHT)-((CELL)*30.0/width));

	glEnd();

	double obj_rwidth=15.0*rwidth;
	double obj_rheight=15.0*rheight;
	double obj_step=(1.0/max)*30.0;
	double obj_x;
	double obj_y;
	int y_step=(*T).get_step(1);
	glBegin(GL_LINES);
	glColor3f(0.0,0.0, 0.0);

 	for(int y=CELL; y<height;y+=CELL){
		index[0]=0; index[1]=y;
		(*T).go(&index);

		obj_y=obj_step*y;
		obj_x=0;
		for(int x=1; x<width;x++){

			valueB=(*T).readpp()*SCALE;
			valueA=(*T).read()*SCALE;
		
			glVertex3f(obj_x-obj_rwidth,valueB,obj_y-obj_rheight);
			obj_x+=obj_step;
			glVertex3f(obj_x-obj_rwidth,valueA,obj_y-obj_rheight);
		}	
	}

	for(int x=CELL; x<width;x+=CELL){
		index[0]=x; index[1]=0;
		(*T).go(&index);

		obj_y=0;
		obj_x=obj_step*x;
		for(int y=1; y<height;y++){

			valueB=(*T).read()*SCALE;
			(*T).it+=y_step;
			valueA=(*T).read()*SCALE;

			glVertex3f(obj_x-obj_rwidth,valueB,obj_y-obj_rheight);
			obj_y+=obj_step;
			glVertex3f(obj_x-obj_rwidth,valueA,obj_y-obj_rheight);
			}
	}
	glEnd();

	glEnable (GL_BLEND); glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glBegin(GL_TRIANGLES);

 	for(int y=CELL; y<height;y+=CELL){
		index[0]=0; index[1]=y;
		(*T).go(&index);

		obj_y=obj_step*y;
		obj_x=0;

		for(int x=CELL; x<width;x+=CELL){

			valueB=(*T).read()*SCALE;
			(*T).it-=y_step*CELL;
			valueC=(*T).read()*SCALE;
			(*T).it+=y_step*CELL+CELL;
			valueA=(*T).read()*SCALE;

			meanValue=(valueA+valueB+valueC)/3;

			if(meanValue>0.0){
				if (meanValue*COLOR_A<1.0) glColor4f(1.0-meanValue*COLOR_C, 1.0-meanValue*COLOR_B, 1.0-meanValue*COLOR_A, 0.5);
				else if (meanValue*COLOR_B<1.0) glColor4f(1.0-meanValue*COLOR_C, 1.0-meanValue*COLOR_B, 0.0, 0.5);
				else if (meanValue*COLOR_C<1.0) glColor4f(1.0-meanValue*COLOR_C, 0.0, 0.0, 0.5);
				else glColor4f(0.0, 0.0, 0.0, 0.5);
			}
			else{
				if (meanValue*COLOR_A>-1.0) glColor4f(1.0+meanValue*COLOR_A, 1.0+meanValue*COLOR_B, 1.0+meanValue*COLOR_C, 0.5);
				else if (meanValue*COLOR_B>-1.0) glColor4f(0.0, 1.0+meanValue*COLOR_B, 1.0+meanValue*COLOR_C, 0.5);
				else if (meanValue*COLOR_C>-1.0) glColor4f(0.0, 0.0, 1.0+meanValue*COLOR_C, 0.5);
				else glColor4f(0.0, 0.0, 0.0, 0.5);
			}

			glVertex3f(obj_x-obj_rwidth,valueB,obj_y-obj_rheight);
			glVertex3f(obj_x-obj_rwidth,valueC,obj_y-obj_step*CELL-obj_rheight);
			obj_x+=obj_step*CELL;
			glVertex3f(obj_x-obj_rwidth,valueA,obj_y-obj_rheight);
		}	
	}

	for(int y=0; y<(height-CELL);y+=CELL){
		index[0]=0; index[1]=y;
		(*T).go(&index);

		obj_y=obj_step*y;
		obj_x=0;

		for(int x=0; x<(width-CELL);x+=CELL){

			valueA=(*T).read()*SCALE;
			(*T).it+=CELL;
			valueC=(*T).read()*SCALE;
			(*T).it+=y_step*CELL;
			valueB=(*T).read()*SCALE;
			(*T).it-=y_step*CELL;

			meanValue=(valueA+valueB+valueC)/3;

			if(meanValue>0.0){
				if (meanValue*COLOR_A<1.0) glColor4f(1.0-meanValue*COLOR_C, 1.0-meanValue*COLOR_B, 1.0-meanValue*COLOR_A, 0.5);
				else if (meanValue*COLOR_B<1.0) glColor4f(1.0-meanValue*COLOR_C, 1.0-meanValue*COLOR_B, 0.0, 0.5);
				else if (meanValue*COLOR_C<1.0) glColor4f(1.0-meanValue*COLOR_C, 0.0, 0.0, 0.5);
				else glColor4f(0.0, 0.0, 0.0, 0.5);
			}
			else {
				if (meanValue*COLOR_A>-1.0) glColor4f(1.0+meanValue*COLOR_A, 1.0+meanValue*COLOR_B, 1.0+meanValue*COLOR_C, 0.5);
				else if (meanValue*COLOR_B>-1.0) glColor4f(0.0, 1.0+meanValue*COLOR_B, 1.0+meanValue*COLOR_C, 0.5);
				else if (meanValue*COLOR_C>-1.0) glColor4f(0.0, 0.0, 1.0+meanValue*COLOR_C, 0.5);
				else glColor4f(0.0, 0.0, 0.0, 0.5);
			}

			glVertex3f(obj_x-obj_rwidth,valueA,obj_y-obj_rheight);
			obj_x+=obj_step*CELL;
			glVertex3f(obj_x-obj_rwidth,valueB,obj_y-obj_rheight+obj_step*CELL);
			glVertex3f(obj_x-obj_rwidth,valueC,obj_y-obj_rheight);
		}
	}
	glEnd();

	glPopMatrix();
	glutSwapBuffers();	
}

void setLIST(LIST * Set){
	list=Set;
}

GLvoid DrawAxis(double X, double Y, double Z, double xtic, double ytic, double ztic){

	glColor3f(0.0, 0.0, 0.0);
	
	if (X>0){
	//Draw X axis
		glVertex3f(-15.0*X,-15.0*Y, -15.0*Z);
		glVertex3f(+15.0*X,-15.0*Y, -15.0*Z);
	//Draw Tick Marks:
		if(xtic>0){
			for (double x=-X+xtic; x<X; x+=xtic){
				glVertex3f(15.0*x,-15.0*Y, -15.0*Z);
				glVertex3f(15.0*x,-15.0*Y-0.5, -15.0*Z);
				glEnd();
				glPushMatrix();
				char text[120];
				sprintf(text, "%f", ((x+X)/X)*15.0/list->x_unit); 	
				glTranslatef(15.0*x,-15.0*Y-0.6, -15.0*Z);
				glRotatef(-45.0,0,0,1.0);
				glScalef(0.005,0.005,0.005);
				for (char *p = &text[0]; *p; p++)
					glutStrokeCharacter(GLUT_STROKE_ROMAN,*p);
				glPopMatrix();
				glBegin(GL_LINES);
			}
		}
	}
	//Draw Y axis
	if (Y>0){
		glVertex3f(-15.0*X,-15.0*Y, -15.0*Z);
		glVertex3f(-15.0*X,+15.0*Y, -15.0*Z);
		if(ytic>0){
			for (double y=-Y+ytic; y<Y; y+=ytic){
				glVertex3f(-15.0*X,15.0*y, -15.0*Z);
				glVertex3f(-15.0*X-0.5,15.0*y, -15.0*Z);
				glEnd();
				glPushMatrix();
				char text[120];
				sprintf(text, "%f", (((y+Y)/Y)*15.0)/list->y_unit); 
				glTranslatef(-15.0*X-0.5-4.0,15.0*y, -15.0*Z);
				glScalef(0.005,0.005,0.005);
				for (char *p = &text[0]; *p; p++)
					glutStrokeCharacter(GLUT_STROKE_ROMAN,*p);
				glPopMatrix();
				glBegin(GL_LINES);
			}
		}
	}
	//Draw Z axis
	if (Z>0){
		glVertex3f(-15.0*X,-15.0*Y, -15.0*Z);
		glVertex3f(-15.0*X,-15.0*Y, +15.0*Z);

		if(ztic>0){
			for (double z=-Z+ztic; z<Z; z+=ztic){
				glVertex3f(-15.0*X,-15.0*Y, 15.0*z);
				glVertex3f(-15.0*X,-15.0*Y-0.5, 15.0*z);
				glEnd();
				glPushMatrix();
				glTranslatef(-15.0*X,-15.0*Y-0.6, 15.0*z);
				glRotatef(-45.0,0.0,0.0,1.0);
				glScalef(0.005,0.005,0.005);
				char text[120];
				sprintf(text, "%f", ((z+Z)/Z*15.0)/list->z_unit); 
				for (char *p = &text[0]; *p; p++)
					glutStrokeCharacter(GLUT_STROKE_ROMAN,*p);
				glPopMatrix();
				glBegin(GL_LINES);
			}
		}
	}
}

GLvoid DrawFlatAxis(double X, double Y, double Z, double xtic, double ytic, double ztic){

	glColor3f(0.0, 0.0, 0.0);
	double draw_x, draw_y, draw_xtic, draw_ytic;
	if(x_proj){
		draw_x=Y;
		draw_y=Z;
		draw_xtic=xtic;
		draw_ytic=ytic;
	}
	if(y_proj){
		draw_x=X;
		draw_y=Z;
		draw_xtic=xtic;
		draw_ytic=ztic;
	}
	if(z_proj){
		draw_x=X;
		draw_y=Y;
		draw_xtic=xtic;
		draw_ytic=ytic;
	}

	//Draw X axis
	glVertex2f(-15.0*draw_x,-15.0*draw_y);
	glVertex2f(+15.0*draw_x,-15.0*draw_y);
	//Draw Tick Marks:
	if(xtic>0){
		for (double x=-draw_x+draw_xtic; x<draw_x; x+=draw_xtic){
			glVertex2f(15.0*x,-15.0*draw_y);
			glVertex2f(15.0*x,-15.0*draw_y-0.5);
			glEnd();
			glRasterPos2f(15.0*x-.2, -15.0*draw_y-1.8);

			char text[120];
			sprintf(text, "%2.2f", ((x+X)/X*15.0)/list->x_unit); 
			for (char *p = &text[0]; *p; p++)
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*p);

//			glutBitmapCharacter(GLUT_BITMAP_8_BY_13 ,text);
			glBegin(GL_LINES);
		}
	}
	//Draw Y axis
	glVertex2f(-15.0*draw_x,-15.0*draw_y);
	glVertex2f(-15.0*draw_x,+15.0*draw_y);
	if(ytic>0){
		for (double y=-draw_y+draw_ytic; y<draw_y; y+=draw_ytic){
			glVertex2f(-15.0*draw_x,15.0*y);
			glVertex2f(-15.0*draw_x-0.5,15.0*y);
			glEnd();
			glRasterPos2f(-15.0*draw_x-.9, 15.0*y-.4);
			char text[120];

			sprintf(text, "%2.2f", ((y+Y)/Y*15.0)/list->y_unit); 
			for (char *p = &text[0]; *p; p++)
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*p);

			glBegin(GL_LINES);
		}
	}
}

void StartOrth(void){
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, Xsize, 0, Ysize);
	glScalef(24, -24, 24);
	glTranslatef(17*Y_length,-17*X_length,0);
	glRotatef(-90.0,0.0,0.0,1.0);
	glMatrixMode(GL_MODELVIEW);
}

void StartPers(void){
	glPushMatrix();
	glTranslatef(x_trans,y_trans,-dist);
	glRotatef(xangle,0.0,1.0,0.0);
	glRotatef(yangle,1.0,0.0,0.0);
}

void EndOrth(void){
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

GLvoid DrawGLSceneList()
{
	ELEMENT *temp_e;
	PAIR *temp_j;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	/* Clear The Screen And The Depth Buffer */

	glLoadIdentity();
	if (mode==ORTH){
		StartOrth();


//		glEnable(GL_ALPHA_TEST);
//		glEnable(GL_LINE_SMOOTH);
//		glLineWidthx(4.5f);

		glBegin(GL_LINES);

		DrawFlatAxis (X_length,Y_length,Z_length,X_length*0.1,Y_length*0.1,Z_length*0.1);

		temp_e=list->start;
		while(temp_e!=NULL){
			temp_j=temp_e->start;
			while(temp_j!=NULL){
				glColor3f(temp_j->R,temp_j->B,temp_j->G);
				double x0,y0,x1,y1;
				if(x_proj){
					x0=temp_e->y*Y_length;
					y0=temp_e->z;
					x1=temp_j->element->y*Y_length;
					y1=temp_j->element->z;
				}
				if(y_proj){
					x0=temp_e->x*X_length;
					y0=temp_e->z;
					x1=temp_j->element->x*X_length;
					y1=temp_j->element->z;
				}
				if(z_proj){
					x0=temp_e->x*X_length;
					y0=temp_e->y*Y_length;
					x1=temp_j->element->x*X_length;
					y1=temp_j->element->y*Y_length;
				}
				if(temp_j->wire) glVertex2f(x0,y0);
				if(temp_j->wire) glVertex2f(x1,y1);
				temp_j=temp_j->next;
			}
			temp_e=temp_e->next;
		}	
		glEnd();
	}
	else if (mode==PERS){
		StartPers();
		glBegin(GL_LINES);
		glLineWidth(10.5f);

		DrawAxis (X_length,Y_length,Z_length,X_length*0.1,Y_length*0.1,Z_length*0.1);

		temp_e=list->start;
		while(temp_e!=NULL){
			temp_j=temp_e->start;
			while(temp_j!=NULL){
				glColor3f(temp_j->R,temp_j->B,temp_j->G);
	
				if(temp_j->wire) glVertex3f(temp_e->x*X_length,temp_e->y*Y_length,temp_e->z);
				if(temp_j->wire) glVertex3f(temp_j->element->x*X_length,temp_j->element->y*Y_length,temp_j->element->z);
				temp_j=temp_j->next;
			}
			temp_e=temp_e->next;
		}	
		glEnd();

		glEnable (GL_BLEND); 
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//	glEnable (GL_LIGHTING);

		glBegin(GL_QUADS);
		temp_e=list->start;
		while(temp_e->next!=NULL){
			temp_j=temp_e->start;
			if(temp_e->quad){
				glColor4f(temp_e->R,temp_e->B,temp_e->G, temp_e->A);
				glVertex3f(temp_e->x*X_length,temp_e->y*Y_length,temp_e->z);
				while(temp_j!=NULL){
					glVertex3f(temp_j->element->x*X_length,temp_j->element->y*Y_length,temp_j->element->z);
					temp_j=temp_j->next;
				}
			}
			temp_e=temp_e->next;
		}	
		glEnd();
	}

	if(mode==ORTH) EndOrth();
	glPopMatrix();

	glutSwapBuffers();	
}

void tick(int value) {
	glutPostRedisplay();
	glutTimerFunc(2, tick, 0);
}

void Start(void){
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	glutInitWindowSize(Xsize,Ysize);
	glutInitWindowPosition(0,0);

	window = glutCreateWindow("Matt's Graphic Library");
	InitGL(Xsize,Ysize);
	
	float sizes, increment;

	glEnable(GL_BLEND);
	glGetFloatv(GL_LINE_WIDTH_RANGE, &sizes);
	glGetFloatv(GL_LINE_WIDTH_GRANULARITY, &increment);

	glutDisplayFunc(DrawGLSceneTensor);
//	glutDisplayFunc(DrawGLSceneList);
	glutTimerFunc(2,tick, 0);
	glutMainLoop();
}
