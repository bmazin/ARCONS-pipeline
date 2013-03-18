/*
To install the plugin (replace "linux" as needed to match the other installed plugins):
NOTE: You might need a -fPIC flag on g++ depending on how tempo2 was compiled...

g++ arcons_plug.C -shared -fPIC -o arcons_linux_plug.t2 -I $TEMPO2/include
cp arcons_linux_plug.t2 $TEMPO2/plugins
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"

//#include "/usr/local/tempo2/include/tempo2.h"

using namespace std;

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
	char parFile[1][MAX_FILELEN];
	char timFile[1][MAX_FILELEN];

	int N_tim = 5000; //Number of photons to process in each batch

	int par_file = 0, a_file = 0, o_file = 0;

	FILE *aFile, *oFile;

	*npsr = 1;

	printf("Arcons Plugin\n");
	printf("Authors: Gualtiero Spiro Jaeger and Michael Johnson\n");
	printf(" --- type ’h’ for help information\n");

	//Read command-line arguments
	for (int i=2;i<argc;i++)
	{
		if (strcmp(argv[i],"-f")==0)
		{
			par_file = 1;
			strcpy(parFile[0],argv[++i]);
		}
		else if (strcmp(argv[i],"-a")==0)
		{
			a_file = 1;
			aFile = fopen(argv[++i],"r");
			if (!aFile)
			{
				fprintf(stderr,"Cannot Open ARCONS File!\n");
				exit(1);
			}
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			o_file = 1;
			oFile = fopen(argv[++i],"w");
			if (!oFile)
			{
				fprintf(stderr,"Cannot Open Output File!\n");
				exit(1);
			}
		}
		else if (strcmp(argv[i],"-h")==0)
		{
			printf("\n TEMPO2 ARCONS plugin\n");
			printf("======================\n");
			printf("\n USAGE: \n\t tempo2 -gr arcons -f par.par -a arcons.toa -o output.txt\n");
			printf("\t -h: this help.\n");
			printf("===============================================\n");
			exit(0);			
		}
	}
    
	if (!par_file)
	{
		printf("No ephemeris (par) file!\n");
		exit(1);
	}

	if (!a_file)
	{
		printf("No ARCONS file!\n");
		exit(1);
	}

	if (!o_file)
	{
		printf("No output file!\n");
		exit(1);
	}

	//Create a dummy timing file
	strcpy(timFile[0],"dummy_timing.tim");
	FILE *f_out;
	f_out = fopen(timFile[0],"w+");
	fprintf(f_out,"FORMAT 1\n");
	fprintf(f_out," dummy 0.0 1.0 0.00000 1\n");
	fclose(f_out);

	//Now read the parfile
	readParfile(psr,parFile,timFile,*npsr);


	//Define our reference line 
	char TZR_tim_str[120];
	sprintf(TZR_tim_str,"FORMAT 1\n TZR %.12Lf %.12Lf 0.0 %s\n",psr->param[param_tzrfrq].val[0], psr->param[param_tzrmjd].val[0], psr->tzrsite);

	char temp_str[120];
	double atoa[N_tim];
	
	int n = 0; 

	int r = 1;

	int tot_photons = 0;

	do
	{
		r = fscanf(aFile,"%s",temp_str);

		if (r == 1)
		{
			atoa[n] = atof(temp_str);
			n++;					
		}

		if (n == N_tim || r != 1)
		{
			f_out = fopen(timFile[0],"w");
			fprintf(f_out,"%s",TZR_tim_str);
			for (int i=0;i<n;i++)
			{
				fprintf(f_out," arcons%d 0.0 %.12f 0.0 pal\n",i,atoa[i]);
			}
			fclose(f_out);

			readTimfile(psr,timFile,*npsr);
			preProcess(psr,*npsr,argc,argv);

			//Form the barycentric arrival times
			formBatsAll(psr,*npsr);

			//Form the residuals (to get the phases)
			formResiduals(psr,*npsr,1);

			for (int i=1;i<n+1;i++)
			{
				double intpart;
				fprintf(oFile,"%.12Lf %f\n",psr[0].obsn[i].bat,modf((double)psr[0].obsn[i].phase,&intpart));
				tot_photons++;
			}

			n=0;
			
			fprintf(stderr,"Finished %d Photons...\n",tot_photons);
		}
	} while (r == 1);

	fclose(aFile);
	fclose(oFile);

	fprintf(stderr,"Finished Calculating BTOAs and Phases for %d Photons\n",tot_photons);
	
	return 0;
}


