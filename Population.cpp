#include "Population.h"
#include "main.h"
#include "Fitness.h"
#include "OutputFormat.h"
#include "Parconst.h"

using namespace std;


// constructors and destructor

Population::Population()
{
}


Population::Population(long int size)
{
    for (long int i = 0; i <= size; i++)
    {
        Individual indiv;
        pop.push_back(indiv);
    }
}


Population::Population(const Population & copy)
    : pop(copy.pop)
{
}


Population::Population(const std::vector<Individual>& vecindiv)
    : pop(vecindiv)
{
}


Population::Population(const ParameterSet& param)
{
    initialize(param);
}



// operator overload

Population & Population::operator=(const Population& copy)
{
    if (this == &copy)
        return(*this);

    pop = copy.pop;
    return(*this);
}


// instance and initialization

void Population::initialize(const ParameterSet& param)
{
    int popsize = param.getpar(INIT_PSIZE)->GetInt();
    //pop.resize(popsize);
    for (long int i = 0; i < popsize; i++)
    {
        Individual indiv(param);
        pop.push_back(indiv);
        //cout << i << endl;
    }
    update();
}


// functions

double fun_sqrt(double x)
{
    return(std::sqrt(x));
}


Population Population::reproduce(long int offspr_number) const
{
    Population offspring;
    vector<double> cumul_fit = cumul_fitness();

    if (offspr_number == 0)
    {
        offspr_number = size();
    }

    offspring.pop.resize(offspr_number);

    for (long int i = 0; i < offspr_number; i++)
    {
        offspring.pop[i] = Individual::mate(
                               this->pick_parent(cumul_fit),
                               this->pick_parent(cumul_fit));
    }
    offspring.update();
    return(offspring);
}


void Population::update(void)
{
    double popvalue = Fitness::GetPopulationValue(*this);
    for (vector<Individual>::iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->update_fitness(popvalue);
    }
}


vector<double> Population::phenotypes() const
{
    vector<double> pheno;
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        pheno.push_back(indiv->get_phenotype());
    }
    return(pheno);
}


double Population::mean_phenotype() const
{
    double sum=0;
    for (vector<Individual>::const_iterator i = pop.begin();
            i != pop.end(); i++)
    {
        sum += i->get_phenotype();
    }
    return(sum/(pop.size()));
}


long int Population::size() const
{
    return(pop.size());
}


vector<double> Population::cumul_fitness() const
{
    vector<double> cum_fit(this->size());
    double cumul = 0;

    vector<double>::iterator cc = cum_fit.begin();
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end();
            indiv++, cc++)
    {
        cumul += indiv->get_fitness();
        *cc = cumul;
    }
    for (vector<double>::iterator i = cum_fit.begin();
            i != cum_fit.end(); i++)
    {
        *i = *i/cumul;
    }
    return(cum_fit);
}


const Individual & Population::pick_parent(const vector<double>& cumfit) const
{
    // return(iterator_search_fit_table(rnum, cumfit));
    return(pop[search_fit_table(Random::randnum(), cumfit)]);
}


long int Population::search_fit_table(double rnum, const vector<double>& cumfit) const
{
    return(sequential_search_fit_table(rnum, cumfit));
}


long int Population::sequential_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    long int i = 0;
    while (cumfit[i] < rnum)
        i++;
    return(i);
}


Individual Population::iterator_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    // Does not work!
    // return(*(std::find_if(cumfit.begin(), cumfit.end(), std::bind2nd(std::less<double>(), rnum))));
    return(Individual());
}


void Population::draw_mutation()
{
    for (unsigned int i = 0; i < pop.size(); i++) {
        pop[i].draw_mutation();
    }
    // population.pop.draw_mutation(population.pop);
}


void Population::make_mutation()
{
    int ind = floor(Random::randnum()*pop.size());
    pop[ind].make_mutation();
}


// output

void Population::write() const
{
    if (!OutputFormat::isInitialized())
    {
        cerr << "Warning: No output!\n";
    }
    write_debug (OutputFormat::GetDebug());
    write_xml   (OutputFormat::GetXml());
    write_simple(OutputFormat::GetSimple());
    write_summary(OutputFormat::GetSummary());
}


void Population::write_debug(ostream & out) const
{
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->write_debug(out);
    }
}


void Population::write_xml(ostream & out) const
{
    out << "xml output: not implemented yet.\n";
}


void Population::write_simple(ostream & out) const
{
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->write_simple(out);
    }
}


void Population::write_summary(ostream & out) const
{
    double sumphen = 0.0, sumphen2 = 0.0;
    double sumfit = 0.0, sumfit2 = 0.0;

    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        sumfit += indiv->get_fitness();
        sumfit2 += indiv->get_fitness()*indiv->get_fitness();

        sumphen += indiv->get_phenotype();
        sumphen2 += indiv->get_phenotype()*indiv->get_phenotype();
    }

    double N = double(pop.size());
    out << "MeanPhen = " << sumphen/N << "\t";
    out << "VarPhen = " << sumphen2/N - (sumphen/N)*(sumphen/N) << "\t";
    out << "MeanFit = " << sumfit/N << "\t";
    out << "VarFit = " << sumfit2/N - (sumfit/N)*(sumfit/N) << "\t";
    out << "FitOpt = " << Fitness::current_optimum() << "\t";
    out << endl;
}

