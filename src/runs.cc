/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <cstdio>
#include <unistd.h>
#include <sys/wait.h>
#include <libgen.h>
#include <cerrno>
#include <cmath>
#include <vector>
#include "../gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

static const char* progname;

typedef struct _process {
	pid_t pid;
	int pipefd[2];
} process;

typedef struct _record {
	double kT;
	double M;
	double M_err;
	double E;
	double E_err;
	double variance; /*of individual lattice point energy*/
	int steps;
} record;

bool record_cmp(record a, record b)
{
	return a.kT < b.kT;
}

pid_t test_pid;
bool is_same_process (const process &other)
{
	return test_pid == other.pid;
}
bool (* is_same_process_gen (const pid_t pid))(const process&)
{
	test_pid = pid;
	return is_same_process;
}

/* Print some documentation. */
void usage()
{
	cout << "Usage: " <<progname<< " [options]\n"\
		"Simulate a 2D ferromagnet using a monte-carlo ising model.\n\n"\
		"  -t time          Run each temperature run for time iterations.\n"\
		"  -s size          Use a size by size lattice.\n"\
		"  -j n             Use n processes (best set to number of CPUs)\n"\
		"  -J J             Set interation parameter to J.\n"\
		"  -H H             Set external Magnetic field to H.\n"\
		"  -a n             Calculate variances over last n iterations.\n"
		"  -o output.dat    Save the data to output.dat.\n"\
		"  -g graph.png     Draw a graph of magnetisation and energy to graph.png.\n"\
		"  -d graph.png     Plot a discrete heat capacity estimate to graph.png\n"\
		"  -f graph.png     Plot a estimate of heat capacity from fluctuations to graph.png\n"\
		"  -1 temp          Lowest temperature to simulate.\n"\
		"  -2 temp          Highest temperature to simulate.\n"\
		"  -n number        Number of temperature steps to take.\n"\
		"  -p               Plot only (reuse data)\n"\
		"  -P               Plot pdf (default png)\n"\
		"  -h               Show this help text.\n";
	exit(1);
}

int main(int argc, char ** argv)
{
	int opt, num_temps=1,parallel=1,status=0;
	bool d=false,f=false,p=false,pdf=false;
	vector<process> children;
	vector<record> records;
	record next_record;
	pid_t tempid;
	const char *output_filename="ising.dat", *heat_filename = NULL, *graph_filename=NULL;
	double kTfrom=DEFAULT_kT, kTto=DEFAULT_kT,kT;
	vector<char *> parameters;
	ofstream output;
	progname = argv[0];
	char command[255];
	command[0] = '\0';
	strncat(command,dirname(argv[0]),255);
	strncat(command,"/ising",255);
	parameters.push_back(command);
	/* Read command line options. */
	while ((opt = getopt(argc,argv,"a:g:d:f:t:s:J:H:o:pP1:2:j:n:h?")) != -1)
	{
		switch(opt)
		{
			case 'a':
				if (atoi(optarg) < 2)
				{
					cout << "Must average over at least 2 iterations!" << endl;
					exit(2);
				}
				parameters.push_back("-a");
				parameters.push_back(optarg);
				break;
			case 's':
				if (atoi(optarg) < 1)
				{
					cout << "Lattice too small!" << endl;
					exit(2);
				}
				parameters.push_back("-s");
				parameters.push_back(optarg);
				break;
			case 'n':
				num_temps = atoi(optarg);
				if (num_temps < 1)
				{
					cout << "Must have at least 1 temperature!" << endl;
					exit(2);
				}
				records.reserve(num_temps);
				break;
			case 'j':
				parallel = atoi(optarg);
				if (parallel < 1)
				{
					cout << "Must have at least 1 subprocess!" << endl;
					exit(2);
				}
				break;
			case '1':
				kTfrom = atof(optarg);
				if (kTfrom < 0)
				{
					cout << "Temperature below 0 is meaningless!" << endl;
					exit(2);
				}
				break;
			case '2':
				kTto = atof(optarg);
				if (kTto < 0)
				{
					cout << "Temperature below 0 is meaningless!" << endl;
					exit(2);
				}
				break;
			case 't':
				if (atoi(optarg) < 1)
				{
					cout << "Must have at least 1 step!" << endl;
					exit(2);
				}
				parameters.push_back("-t");
				parameters.push_back(optarg);
				break;
			case 'J':
				parameters.push_back("-J");
				parameters.push_back(optarg);
				break;
			case 'H':
				parameters.push_back("-H");
				parameters.push_back(optarg);
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 'g':
				graph_filename = optarg;
				break;
			case 'd':
				heat_filename = optarg;
				d = true;
				break;
			case 'f':
				heat_filename = optarg;
				f = true;
				break;
			case 'p':
				p = true;
				break;
			case 'P':
				pdf = true;
				break;
			case 'h':
			case '?':
			default:
				usage();
		}
	}

	if (p)
		goto plot;
	
	/* We need from and to in the right order or later a loop will go on forever */
	if (kTfrom>kTto)
	{
		double kTtmp = kTto;
		kTto = kTfrom;
		kTfrom = kTtmp;
		cerr << "Swapping temperatures 1 and 2." << endl;
	}

	output.open(output_filename);

	/* Do a run at each temperature. We detect kTfrom==kTto later*/
	for (kT=kTfrom; kT <= kTto; kT += (kTto-kTfrom)/num_temps)
	{
		if (children.size() >= parallel)
		{
			tempid = wait(&status);
			if (tempid == -1)
			{
				perror("Could not wait for child");
				exit(1);
			}
			if (status != 0)
			{
				if (WIFEXITED(status))
					cerr << "A child exited with status " << WEXITSTATUS(status)<<endl;
				else if (WIFSIGNALED(status))
					cerr << "A child was killed by signal " << WTERMSIG(status)<<endl;
				else if (WIFSTOPPED(status))
					cerr << "A child was stopped  by signal " << WSTOPSIG(status)\
					<< ".\n Don't do that." << endl;
				else if (WIFCONTINUED(status))
					cerr << "A child was continued (?)" <<endl;
				else
					cerr << "A child went missing with status " << status << endl;
				exit(1);
			}
			vector<process>::iterator this_process = find_if(children.begin(),
			                        children.end(),is_same_process_gen(tempid));
			fscanf(fdopen(this_process->pipefd[0],"r"),
			       "%lf %lf %lf %lf %lf %lf %d\n",
			       &next_record.kT, &next_record.M, &next_record.M_err,
			       &next_record.E, &next_record.E_err, &next_record.variance,
			       &next_record.steps);
			records.push_back(next_record);
			close(this_process->pipefd[0]);
			children.erase(this_process);
		}

		process new_process;
		pipe(new_process.pipefd);
		if((tempid=fork())==0)
		{

			if(dup2(new_process.pipefd[1],STDOUT_FILENO)==-1)
			{
				perror("Could not reassign stdout to pipe");
				exit(1);
			}
			close(new_process.pipefd[0]);
			close(new_process.pipefd[1]);
			char tempstr[31];
			sprintf(tempstr,"%.15E",kT);
			parameters.push_back("-T");
			parameters.push_back(tempstr);
			parameters.push_back(NULL);
			execv(command,&parameters[0]);
			perror("Could not run subprocess");
			exit(1);
		}
		if (tempid == -1)
		{
			perror("Could not fork");
			exit(1);
		}
		close(new_process.pipefd[1]);
		new_process.pid = tempid;
		children.push_back(new_process);
		
		/* if kTto==kTfrom the kT increment is zero, and we therefore must only
		 * do one run, so break in this case. */
		if (kTto == kTfrom)
			break;
	}

	for(vector<process>::iterator it = children.begin(); it != children.end(); ++it)
	{
		if(waitpid(it->pid,NULL,0)==-1)
		{
			perror("Could not collect a child (ignoring)");
		}
		fscanf(fdopen(it->pipefd[0],"r"),
		       "%lf %lf %lf %lf %lf %lf %d\n",
		       &next_record.kT, &next_record.M, &next_record.M_err,
		       &next_record.E, &next_record.E_err, &next_record.variance,
		       &next_record.steps);
		records.push_back(next_record);
		close(it->pipefd[0]);
		close(it->pipefd[1]);
	}

	/* We can get out of order when running with parallel > 1 */
	if (parallel > 1)
		sort(records.begin(),records.end(),record_cmp);
	/* Print out output */
	for(vector<record>::iterator it = records.begin(); it != records.end();++it)
	{
		output << it->kT << ' ' << fabs(it->M) << ' ' << it->M_err \
			<< ' ' << it->E << ' ' << it->E_err \
			<< ' ' << (it == records.begin() ? 0.0 : \
				           (it->E-(it-1)->E)/(it->kT-(it-1)->kT)) \
			<< ' ' << it->variance/(it->kT*it->kT)\
			<< ' ' << it->steps << endl;
	}
	
	/* Close the output to ensure it is flushed for gnuplot to read. */
	output.close();

plot:
	/* Tell gnuplot to draw us a graph.*/
	if (graph_filename != NULL)
	{
		Gnuplot gp;
		gp << \
			"set term " << (pdf?"pdf\n":"png size 1024,768\n") << \
			"set output '" << graph_filename << "'\n"\
			"set xlabel 'kT/J'\n"\
			"set ylabel 'Magnetization'\n"\
			"set yrange [0:1]\n"\
			"set xtics rotate by -45 add ('Tc(Onsager)' 2.269185314)\n"\
			"set y2label 'Energy per lattice point/J'\n"\
			"set y2range [-2:0]\n"\
			"set ytics nomirror\n"\
			"set y2tics\n"\
			"plot '"<<output_filename \
			<< "' u 1:2:3 w yerrorbars t 'Magnetization'"\
			",'"<<output_filename << \
			"' u 1:4:5 w yerrorbars t 'Energy' axes x1y2\n";
	}
	if(d || f)
	{
		Gnuplot gp;
		gp << \
			"set term " << (pdf?"pdf\n":"png size 1024,768\n") << \
			"set output '" << heat_filename << "'\n"\
			"set xlabel 'kT/J'\n"\
			"set ylabel 'Energy per lattice point/J'\n"\
			"set xtics rotate by -45 add ('Tc(Onsager)' 2.269185314)\n"\
			"set y2label 'Heat Capacity per lattice point/k'\n"\
			"set ytics nomirror\n"\
			"set yrange [-2:0]\n"\
			"set y2range [0:]\n"\
			"set y2tics\n"\
			"plot '"<<output_filename << "' u 1:4:5 w yerrorbars t 'Energy'";
		if(d)
			gp << ",'"<<output_filename << \
			"' u 1:6 t 'Heat capacity (discrete)' axes x1y2";
		if(f)
			gp << ",'"<<output_filename << \
			"' u 1:7 t 'Heat capacity (fluctuations)' axes x1y2";
		gp << "\n";
	}

	return 0;
}