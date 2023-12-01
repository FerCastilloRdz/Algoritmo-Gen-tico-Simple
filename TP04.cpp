    /*
    Algoritmo generico simple
    */
    //Librerias
    #include<cstdio>
    #include<cstdlib>
    #include<cmath>
    #include<ctime>

    typedef struct{
                unsigned char *Crom; //Cromosoma(cadena de bits) [1 byte or bit del cromosoma]
                double *Vreal; //Vector de valores reales
                int *Vent; //Vector de valor entero
                double VObj; //Valor Objetivo
                double Vfit;//valor de adaptacion
    }INDIVIDUO;

    //Objetos a clase para el algrodimo genetico
    class GA
    {
    private:
        INDIVIDUO *Pob;//Poblacion
        INDIVIDUO *NewPob;//Poblacion
        unsigned int PobSize;//Tamaño de la poblacion
        unsigned int NGens; //Numero de genes
        const unsigned int *BitxGens;//Numero de bits de cada gen
        unsigned int CromSize; //Tamaño total del coromosoma
        const double *LSub;//Limite superiores para decodificar real
        const double *LInf;//Limite inferior para decodificar real
        double SumaObj;//Suma de todos los valores objetivo
        double SumaFit;//
        double PromedioObj;//Promedio de los valores objetivo
        unsigned int IdBest; //indice al mejor individuo de toda la poblacion
        unsigned int IdWorst;//indice del peor individuo de toda la poblacion
        unsigned int *seleccinado;//
        unsigned int maxGeneracion;

    public:
        GA(unsigned int numeroDeGEneraciones,unsigned int SizePob, unsigned int NumDeGenes, const unsigned *NBitXGen, const double *LimitSup, const double *LimitInf)
        {
            maxGeneracion=numeroDeGEneraciones;
            PobSize =  SizePob;
            NGens = NumDeGenes;
            BitxGens = NBitXGen;
            LSub = LimitSup;
            LInf = LimitInf;
            CromSize = 0;
            SumaObj = 0;
            PromedioObj = 0;
            IdBest = 0;
            IdWorst = 0;
            for (int i = 0; i < NGens; ++i)
            {
                CromSize += BitxGens[i];
            }

            //Generar la memoria para una poblacion de N individos
            Pob = new INDIVIDUO[PobSize];
            //Generar memoria para cada individuo
            for(int i = 0; i < PobSize; i++ )
            {
                Pob[i].Crom = new unsigned char [CromSize];
                Pob[i].Vreal = new double [NGens];
                Pob[i].Vent = new int [NGens];
            }

            seleccinado=new unsigned int[PobSize];

            NewPob = new INDIVIDUO[PobSize];
            for(int i = 0; i < PobSize; i++ )
            {
                NewPob[i].Crom = new unsigned char [CromSize];
                NewPob[i].Vreal = new double [NGens];
                NewPob[i].Vent = new int [NGens];
            }
            //Inicializar el cromosoma
            for (int k = 0; k < PobSize; ++k)//Para toda la poblacion
            {
                for (int i = 0; i < CromSize; ++i)//Para todo el cromosoma del poblacor K-esimo inidividuo
                {
                    Pob[k].Crom[i] = rand()%2;
                }
            }
        }

        void ImprimeInd(unsigned int id)
        {
            unsigned int acumulador,gen;
            //Impribe el cromosoma del individuo id
            printf("\n%i: ", id);

            gen=0;
            acumulador = BitxGens[gen];

            for(int k = 0; k < CromSize; k++)
            {
                if(k ==acumulador)
                {
                    printf(" : ");
                    gen++;
                    acumulador += BitxGens[gen];
                }

                 printf("%i",Pob[id].Crom[k]);
            }
            printf(" [");
            for (int k = 0; k < NGens; ++k)
            {
                printf(" %3i ",Pob[id].Vent[k]);
            }
            printf("] ");
            printf(" [");
            for (int i = 0; i < NGens; ++i)
            {
                printf("%2.3f, ",Pob[id].Vreal[i]);
            }
            printf("] ");
            printf(" [");
            printf(" Vobj:%2.3f ",Pob[id].VObj);
            printf("] ");
            printf(" [");
            printf(" FIT:%2.3f ",Pob[id].Vfit);
            printf("] ");

        }
        void ImpirmeResumen()
        {
            printf("\nBest = %i, Worst = %i", IdBest, IdWorst);
            printf("\nSuma obj:%2.3f, ",SumaObj);
            printf("Promedio obj:%2.3f, ",PromedioObj);
        }
        void ImprimeBest()
        {
            ImprimeInd(IdBest);
        }
        void ImpirmePoblacion()
        {
            for (int i = 0; i < PobSize; ++i)//Para toda la poblacion
            {
                ImprimeInd(i);
            }
        }

        void DecodeEntero()
        {
            unsigned int acumulador, pos;
            for (int k = 0; k < PobSize; ++k)//Para toda la poblacon
            {
                int j = 0;
                acumulador = 0;

                for (int i = 0; i < NGens; ++i)//Para cada Gen del K-esimo del individuo
                {
                    Pob[k].Vent[i] = 0;
                    acumulador += BitxGens[i];
                    pos = 0;

                    for ( ; j < acumulador; ++j)//Para cada bit de cada gen
                    {
                        Pob[k].Vent[i] += Pob[k].Crom[j] * pow(2,pos);
                        pos++;
                    }

                }
            }
        }

        void DecodeReal()//Decodificar un valor real
        {
            double porcentaje, rango;
            DecodeEntero();
            for (int k = 0; k < PobSize; ++k)
            {
                for (int i = 0; i < NGens; ++i)
                {
                    rango = LSub[i]-LInf[i];
                    porcentaje = Pob[k].Vent[i]/(pow(2,BitxGens[i])-1);
                    Pob[k].Vreal[i] = (porcentaje * rango) + LInf[i];
                }
            }
        }

        double FuncionObjetivo(unsigned int i)
        {
        double fit=0;
        int k;
        //Pob[i].Vreal[0]=0;
        for(k=0;k<30;k++)
        {
        fit+=(pow(Pob[i].Vreal[0],2)-10*cos(6.2831*Pob[i].Vreal[0])+10);
        }
    //    printf("\n El resultado es %f",fit);
    //    exit(0);
        return (fit);
        }
        void Obj_a_fit()
        {
            unsigned int k;
            double rango,porcentaje;
            rango=Pob[IdBest].VObj-Pob[IdWorst].VObj;
           if(fabs(rango)<0.000001)
           {SumaFit=PobSize*100;
            for(k=0;k<PobSize;k++)
            {
                Pob[k].Vfit = 100;
            }
           }
           else
           {SumaFit=0;
            for(k=0;k<PobSize;k++)
            {
                porcentaje = (Pob[k].VObj-Pob[IdWorst].VObj)/rango;
                Pob[k].Vfit = 100*porcentaje;
                SumaFit+=Pob[k].Vfit;
            }
           }
        }


       void NextGeneracion()
       {
           INDIVIDUO *aux;
           aux=Pob;
           Pob=NewPob;
           NewPob=aux;
       }

        void seleccion()
        {
         double *Ruleta = new double[PobSize];
         double suma=0,suma1=0,pelota;
         //crear la ruleta
        for(int k=0;k<PobSize;k++)
         {
             Ruleta[k]=Pob[k].Vfit/SumaFit;

             suma1=Ruleta[k];

             //printf("\nRuleta [%i] = %g %c",k,suma1*100,37);

             suma+=Ruleta[k];
         }
         //printf("\nla suma de la ruleta VFIT es = %g %c",suma*100,37);

         for(int i=0;i<PobSize;i++)
        {
          pelota=(double)rand()/RAND_MAX;
          suma=0;
          for(int k=0;k<PobSize;k++)
           {
              suma+=Ruleta[k];
              if(suma>pelota)
              {
                  seleccinado[i]=k;
                  break;
              }
           }
           //printf("\nSeleccionado[%i]: %i",i,seleccinado[i]);
         }

        }

        void cruza(double Pc)
        {
            int PdC,b;
            double r;
          for(int k=0;k<PobSize;k+=2)
         {  r=(double)rand()/RAND_MAX;
            if(r<Pc)
            {//se aplica la cruza
                PdC=rand()%CromSize;
                //printf("\n %i: punto de cruza= %i",k,PdC);
                for(int i=0;i<=PdC;i++)
                {
                    NewPob[k].Crom[i]=Pob[seleccinado[k]].Crom[i];
                    NewPob[k+1].Crom[i]=Pob[seleccinado[k+1]].Crom[i];
                }
                for(int i=PdC+1;i<CromSize;i++)
                {
                    NewPob[k+1].Crom[i]=Pob[seleccinado[k]].Crom[i];
                    NewPob[k].Crom[i]=Pob[seleccinado[k+1]].Crom[i];
                }

            }
            else
            {//no se aplica la cruza
                for(int i=0;i<CromSize;i++)
                {
                    NewPob[k].Crom[i]=Pob[seleccinado[k]].Crom[i];
                    NewPob[k+1].Crom[i]=Pob[seleccinado[k+1]].Crom[i];
                }
            }
         }

        }

        void Muta(double PM)
        {
            int PdC,b;
            double r;
            for(int k=0;k<PobSize;k++)
            {
                for(int i=0;i<CromSize;i++)
                {  r=(double)rand()/RAND_MAX;
                    if(r<PM)
                     {//se aplica la mutacion
                        if(NewPob[k].Crom[i])
                            NewPob[k].Crom[i]=0;
                        else
                            NewPob[k].Crom[i]=1;
                     }

                }
            }
        }
        void Elitismo()
        {
             for(int i=0;i<CromSize;i++)
                {
                    NewPob[0].Crom[i]=Pob[IdBest].Crom[i];
                }
        }
        void EvaluarPob()
        {
            double aux, BestObj,WorstObj;
            unsigned int k;
            IdBest = 0;
            IdWorst = 0;
            SumaObj = 0;
            /*Evaluar el primer individuo o convertir el
               valor objetivo en el mejor y el peor*/
            aux = FuncionObjetivo(0);
            Pob[0].VObj = aux;
            BestObj = aux;
            WorstObj = aux;
            SumaObj += aux;

            for (k = 1; k < PobSize; ++k)//Para toda la poblacion
            {
                aux = FuncionObjetivo(k);
                Pob[k].VObj = aux;
                SumaObj += aux;
                if(aux<BestObj)
                {
                    BestObj = aux;
                    IdBest = k;
                }

                if (aux>WorstObj)
                {
                    WorstObj = aux;
                    IdWorst = k;
                }
            }
            PromedioObj = SumaObj/PobSize;
        }
        double promedio(unsigned int Id)
    {
        double merjoPosi=0;
        merjoPosi = Pob[IdBest].VObj;
        return merjoPosi;
    }
        ~GA()
        {
            //Liberar la memoria para cada individuo
            for(int i = 0; i < PobSize; i++ )
            {
                delete Pob[i].Crom;
                delete Pob[i].Vreal;
                delete Pob[i].Vent;
            }

            //Liberar la memoria para la poblacion
            delete[] Pob;
        };


    };//Fin de la clase GA


    //Parametros del Algoritmos Genetico
    const double numeroDeGEneraciones= 200;
    const unsigned int SizePoblacion = 40;
    const unsigned int NumeroDeGenes = 30; //0-511
    const unsigned int NumerosDeBitsxGen[NumeroDeGenes] = {14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14};
    const double LimiteSuperior[NumeroDeGenes] = { 5.12, 5.12, 5.12, 5.12, 5.12, 5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12}; //De -1 a 1
    const double LimiteInferior[NumeroDeGenes] = {-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};

    int main()
    {
        float promedioDatos=0,sumaDatos=0,desvi=0,sumatoria=0;
        for(int i=1;i<=100;i++)
        {
        srand(time(0));
        GA ga(numeroDeGEneraciones,SizePoblacion, NumeroDeGenes, NumerosDeBitsxGen, LimiteSuperior, LimiteInferior);
        printf("\nSe genero el GA...");
        int t=1;
        ga.DecodeReal();
        ga.EvaluarPob();
        ga.Obj_a_fit();
        printf("\nGeneracion: 0");
        //ga.ImpirmePoblacion();
        //ga.ImpirmeResumen();
        while(t<=numeroDeGEneraciones)
        {
            ga.seleccion();
            ga.cruza(1.0);
            ga.Muta(0.1);
            ga.Elitismo();
            ga.NextGeneracion();
            ga.DecodeReal();
            ga.EvaluarPob();
            ga.Obj_a_fit();
            //printf("\nGeneracion: %i",t);
            //ga.ImpirmePoblacion();
            //ga.ImpirmeResumen();
            t++;
        }
        printf("\nGeneracion: %i",t);
        printf("\nResultado de la ejecucion: %i",i);
        //ga.ImpirmePoblacion();
        ga.ImprimeBest();
        ga.ImpirmeResumen();
        printf("\nmejor valor:%0.2f   \n",ga.promedio(i));
        sumaDatos += ga.promedio(i);
        promedioDatos=sumaDatos/100;
        sumatoria+=pow((ga.promedio(i)-promedioDatos),2);
        desvi=sqrt(sumatoria/100);
        }
        //promedioDatos=sumaDatos/100;
        printf("\n el promedio de las 100 ejecuciones es:  %.3f  ",promedioDatos);
        printf("\n la desviacion estnadar de las 100 ejecuciones es:  %.3f  ",desvi);
        return 0;
    }//FIN DEL MAIN
