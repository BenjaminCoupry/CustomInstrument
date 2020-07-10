using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;
using System.IO;

namespace CustomInstrument
{
    class Program
    {
        static void Main(string[] args)
        {
            Random r = new Random();
            Enveloppe env = new Enveloppe(ProfilEnveloppe.Cos(0, 1, 20000), new ProfilEnveloppe(new List<ProfilEnveloppe>(2) { ProfilEnveloppe.Cos(1, 0.5, 6000), ProfilEnveloppe.Cos(0.5, 1, 9000) }), ProfilEnveloppe.Cos(1,0,30000));
            ProfilHarmoniques ph = ProfilHarmoniques.Alea(ref r, 3);
            Instrument flutador = new Instrument(env, ph, 220,0.8);
            double[] do_ = flutador.Jouer(4, 0, 44100,0);
            string path = "D:/lab/musik/clair.txt";
            //DoubleToFile(k, path+".raw");
            //VersWav(44100, path + "test.wav", do_);
            Orchestre orch= new Orchestre(new Instrument[] { flutador });
            Musique mc = new Musique(44100, orch, Note.McFromFile(path));
            VersWav(44100, path+".musique.wav", mc.Jouer());

        }
        //http://soundfile.sapp.org/doc/WaveFormat/
        static void VersWav(int fe,string Path, double[] inputs)
        {
            uint numsamples = (uint)inputs.GetLength(0);
            ushort numchannels = 1;
            ushort samplelength = 2; // in bytes
            uint samplerate = (uint)fe;

            FileStream f = new FileStream(Path, FileMode.Create);
            BinaryWriter wr = new BinaryWriter(f);
            wr.Write(System.Text.Encoding.ASCII.GetBytes("RIFF"));
            wr.Write((int)(36 + numsamples * numchannels * samplelength));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("WAVE"));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("fmt "));
            wr.Write(16);
            wr.Write((ushort)1);
            wr.Write((ushort)numchannels);
            wr.Write((int)samplerate);
            wr.Write((int)(samplerate * samplelength * numchannels));
            wr.Write((ushort)(samplelength * numchannels));
            wr.Write((ushort)(8 * samplelength));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("data"));
            wr.Write((int)(numsamples * samplelength));
            for (int i = 0; i < numsamples; i++)
            {
                wr.Write(BitConverter.GetBytes((short)(inputs[i]*(double)short.MaxValue)));
            }
        }

        public static void Etendre(ref List<byte> ls, byte[] add)
        {
            for(int i=0;i<add.GetLength(0);i++)
            {
                ls.Add(add[i]);
            }
        }
        public static byte[] StringToByteArray(string hexString)
        {
            if (hexString == null)
                throw new ArgumentNullException("hexString");
            if (hexString.Length % 2 != 0)
                throw new ArgumentException("hexString must have an even length", "hexString");
            var bytes = new byte[hexString.Length / 2];
            for (int i = 0; i < bytes.Length; i++)
            {
                string currentHex = hexString.Substring(i * 2, 2);
                bytes[i] = Convert.ToByte(currentHex, 16);
            }
            return bytes;
        }
        static byte[] DoubleToByte(double[] input)
        {
            byte[] retour = new byte[sizeof(double) * input.GetLength(0)];
            int k = 0;
            for(int i=0; i< input.GetLength(0); i++)
            {
                byte[] val = BitConverter.GetBytes(input[i]);
                for (int j = 0; j < sizeof(double); j++)
                {
                    retour[k] = val[j];
                    k++;
                }
            }
            return retour;
        }
        static void ByteToFile(byte[] input, string Path)
        {
            System.IO.File.WriteAllBytes(Path,input);
        }
        static void DoubleToFile(double[] input, string Path)
        {
            System.IO.File.WriteAllBytes(Path, DoubleToByte(input));
        }
        static double[] TxtRead(string Path)
        {
            string[] lignes = System.IO.File.ReadAllLines(Path);
            double[] retour = new double[lignes.GetLength(0)];
            for(int i=0;i<lignes.GetLength(0);i++)
            {
                retour[i] = Convert.ToDouble(lignes[i]);
            }
            return retour;
        }

        class ProfilHarmoniques
        {
            private int nbHarmoniques;
            private double[] intensites;

            public ProfilHarmoniques(double[] intensites)
            {
                this.intensites = intensites;
                nbHarmoniques = intensites.GetLength(0);
            }

            public bool HarmoniquesValides()
            {
                double somme = 0;
                for(int i=0;i<nbHarmoniques;i++)
                {
                    somme += Math.Abs(intensites[i]);
                }
                return somme <= 1;
            }

            public double GetIntensite(int n)
            {
                return intensites[n];
            }
            public int GetNbHarmoniques()
            {
                return nbHarmoniques;
            }
            public static ProfilHarmoniques Alea(ref Random r, int harmMax)
            {
                int nb = r.Next(1, harmMax + 1);
                double[] intensites = new double[nb];
                for(int i=0;i<nb;i++)
                {
                    intensites[i] = r.NextDouble();
                }
                double sum = intensites.Sum();
                for (int i = 0; i < nb; i++)
                {
                    intensites[i] /= sum;
                }
                return new ProfilHarmoniques(intensites);
            }
        }
        class ProfilEnveloppe
        {
            private int nbEchant;
            private double[] intensites;

            public ProfilEnveloppe(double[] intensites)
            {
                this.intensites = intensites;
                nbEchant = intensites.GetLength(0);
            }
            public ProfilEnveloppe(List<ProfilEnveloppe> profils)
            {
                nbEchant = profils.Sum(e => e.nbEchant);
                intensites = new double[nbEchant];
                int k = 0;
                for(int i=0;i<profils.Count;i++)
                {
                    ProfilEnveloppe pr = profils.ElementAt(i);
                    for(int j=0;j<pr.nbEchant;j++)
                    {
                        intensites[k] = pr.GetIntensite(j);
                        k++;
                    }
                }
            }

            public double GetIntensite(int n)
            {
                return intensites[n];
            }
            public int GetNbEchants()
            {
                return nbEchant;
            }
            public static ProfilEnveloppe Lineaire(double i0, double ifin, int n)
            {
                double[] ret = new double[n];
                for(int i=0;i<n;i++)
                {
                    double k = (double)i / (double)n;
                    ret[i] = i0 + k * (ifin - i0);
                }
                return new ProfilEnveloppe(ret);
            }
            public static ProfilEnveloppe Tremblante(double i0, double ifin,int nbTremblements,double tassement, int n)
            {
                double[] ret = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double k = (double)i / (double)n;
                    ret[i] = (i0 + k * (ifin - i0))*(Math.Cos(2.0*Math.PI*nbTremblements*k)+1+tassement)/(2+tassement);
                }
                return new ProfilEnveloppe(ret);
            }
            public static ProfilEnveloppe Puissance(double i0, double ifin,double beta, int n)
            {
                double[] ret = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double k = Math.Pow((double)i / (double)n,beta);
                    ret[i] = i0 + k * (ifin - i0);
                }
                return new ProfilEnveloppe(ret);
            }
            public static ProfilEnveloppe Alea(double i0, double ifin, int nbAlea, ref Random r, int n)
            {
                double[] retour = new double[n];
                double[] rand = new double[nbAlea + 2];
                rand[0] = i0;
                rand[nbAlea + 1] = ifin;
                for(int i=1;i<nbAlea+1;i++)
                {
                    rand[i] = r.NextDouble();
                }

                int i_ = -1;
                int p = (int)Math.Ceiling((float)n / (float)(nbAlea + 1));
                for (int i=0;i<n;i++)
                {
                    if (i % p == 0)
                    {
                        i_++;
                    }
                    double k = (double)(i % p) / (double)p;
                    double u0 = rand[i_];
                    double u1 = rand[i_ + 1];
                    retour[i] = u0 + 0.5*(1.0-Math.Cos(Math.PI*k)) * (u1 - u0);
                    
                }
                return new ProfilEnveloppe(retour);
            }
            public static ProfilEnveloppe Exp(double i0, double ifin, int n)
            {
                double lambda = Math.Log(ifin / i0);
                double[] ret = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double k = (double)i / (double)n;
                    ret[i] = i0*Math.Exp(lambda*k);
                }
                return new ProfilEnveloppe(ret);
            }
            public static ProfilEnveloppe Cos(double i0, double ifin, int n)
            {
                double[] ret = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double k = (double)i / (double)n;
                    ret[i] = i0 + 0.5*(1.0-Math.Cos(Math.PI*k)) * (ifin - i0);
                }
                return new ProfilEnveloppe(ret);
            }
        }
        class Enveloppe
        {
            private ProfilEnveloppe attaque;
            private ProfilEnveloppe courant;
            private ProfilEnveloppe fin;

            public Enveloppe(ProfilEnveloppe attaque, ProfilEnveloppe courant, ProfilEnveloppe fin)
            {
                this.attaque = attaque;
                this.courant = courant;
                this.fin = fin;
            }
            public double[] GetEnveloppe(int nbEchant)
            {
                int repetCourant = Math.Max(0, (nbEchant - (attaque.GetNbEchants() + fin.GetNbEchants())) / courant.GetNbEchants());
                int echanttotal = attaque.GetNbEchants() + fin.GetNbEchants() + repetCourant * courant.GetNbEchants();
                double[] retour = new double[echanttotal];
                int k = 0;
                for(int i=0;i< attaque.GetNbEchants(); i++)
                {
                    retour[k] = attaque.GetIntensite(i);
                    k++;
                }
                for (int n = 0; n < repetCourant; n++)
                {
                    for (int i = 0; i < courant.GetNbEchants(); i++)
                    {
                        retour[k] = courant.GetIntensite(i);
                        k++;
                    }
                }
                for (int i = 0; i < fin.GetNbEchants(); i++)
                {
                    retour[k] = fin.GetIntensite(i);
                    k++;
                }
                return retour;
            }
            public int GetNbEchantsAttaque()
            {
                return attaque.GetNbEchants();
            }
            public int GetNbEchantsFin()
            {
                return fin.GetNbEchants();
            }
        }
        class Instrument
        {
            private Enveloppe enveloppe;
            private ProfilHarmoniques harmoniques;
            private double f0;
            private double force;

            public Instrument(Enveloppe enveloppe, ProfilHarmoniques harmoniques, double f0, double force)
            {
                this.enveloppe = enveloppe;
                this.harmoniques = harmoniques;
                this.f0 = f0;
                this.force = force;
            }

            public void Accorder(double f_fond)
            {
                f0 = f_fond;
            }
            public double[] Jouer (double T, int note, double fe, int Anticipation)
            {
                int nbEchant = (int)(fe * T);
                nbEchant += Anticipation;
                double fi = f0 * Math.Pow(2, (double)note / 12.0);
                double[] env = enveloppe.GetEnveloppe(nbEchant);
                int n = env.GetLength(0);
                double[] retour = new double[n];
                for(int i=0;i<n;i++)
                {
                    double t = (double)i / fe;
                    double k = 0;
                    for(int j=0;j<harmoniques.GetNbHarmoniques();j++)
                    {
                        k += harmoniques.GetIntensite(j)*Math.Sin(2.0 * Math.PI * (fi * (j + 1)) * t)*env[i];
                    }
                    retour[i] = k*force;
                }
                return retour;
            }
            public Enveloppe GetEnveloppe()
            {
                return enveloppe;
            }
        }
        class Note
        {
            private double T0;
            private double Duree;
            private List<int> IndInstruments;
            public int note;
            private bool Anticiper;
            private bool Complete;
            private double Nuance;

            public Note(double t0, double duree, List<int> indInstruments, int note, bool anticiper, bool complete, double nuance)
            {
                T0 = t0;
                Duree = duree;
                IndInstruments = indInstruments;
                this.note = note;
                Anticiper = anticiper;
                Complete = complete;
                Nuance = nuance;
            }

            public double[] Jouer (Orchestre or, double fe)
            {
                List<double[]> jeux = new List<double[]>();
                int attaqueMax = 0;
                int finMax = 0;
                if (Anticiper || Complete)
                {
                    for (int i = 0; i < IndInstruments.Count; i++)
                    {
                        Instrument Ins = or.GetInstrument(IndInstruments.ElementAt(i));
                        if (Anticiper && Ins.GetEnveloppe().GetNbEchantsAttaque() > attaqueMax)
                        {
                            attaqueMax = Ins.GetEnveloppe().GetNbEchantsAttaque();
                        }
                        if (Complete && Ins.GetEnveloppe().GetNbEchantsFin() > finMax)
                        {
                            finMax = Ins.GetEnveloppe().GetNbEchantsFin();
                        }
                    }
                }
                for (int i = 0; i < IndInstruments.Count; i++)
                {
                    Instrument Ins = or.GetInstrument(IndInstruments.ElementAt(i));
                    jeux.Add(Ins.Jouer(Duree+(double)finMax/fe, note, fe, attaqueMax));
                }
                T0 -= (double)attaqueMax / fe;
                int lgMax = jeux.Max((f) => f.GetLength(0));
                double[] retour = new double[lgMax];
                for(int i=0;i<lgMax;i++)
                {
                    double somme = 0;
                    for (int j = 0; j < IndInstruments.Count; j++)
                    {
                        if(jeux.ElementAt(j).GetLength(0)>i)
                        {
                            somme += jeux.ElementAt(j)[i];
                        }
                    }
                    retour[i] = somme*Nuance;
                }
                return retour;
            }
            public double GetT0()
            {
                return T0;
            }
            public double GetDuree()
            {
                return Duree;
            }
            public static List<Note> McFromFile(string Path)
            {
                List<Note> retour = new List<Note>();
                Dictionary<string, double> Aliases = new Dictionary<string, double>();
                string[] lignes = System.IO.File.ReadAllLines(Path);
                double longTemps = 0;
                double temps = 0;
                for (int i=0;i<lignes.GetLength(0);i++)
                {
                    string ligne = lignes[i];
                    if(ligne.Split()[0]=="TEMPS")
                    {
                        longTemps = Convert.ToDouble(ligne.Split()[1]);
                    }
                    else if(ligne.Split()[0] == "BPM")
                    {
                        longTemps = 60.0/Convert.ToDouble(ligne.Split()[1]);
                    }
                    else if (ligne.Split()[0] == "SKIP")
                    {
                        string st = ligne.Split()[1];
                        double delta;
                        if(double.TryParse(st,out delta))
                        {
                            temps += longTemps * delta;
                        }
                        else
                        {
                            temps += Aliases[st] * longTemps;
                        }

                    }
                    else if (ligne.Split()[0] == "ALIAS")
                    {
                        string cle = ligne.Split()[1];
                        double valeur = Convert.ToDouble(ligne.Split()[2]);
                        if(Aliases.ContainsKey(cle))
                        {
                            Aliases[cle] = valeur;
                        }
                        else
                        {
                            Aliases.Add(cle, valeur);
                        }
                    }
                    else if(ligne =="" || (ligne.Length>=2 && ligne.Substring(0,2)=="//"))
                    {

                    }
                    else if(ligne == ".")
                    {
                        temps += longTemps;
                    }
                    else
                    {
                        string[] infNotes = ligne.Split(new char[] { '+', ' '});
                        foreach(string st in infNotes)
                        {
                            if (st != "")
                            {
                                retour.Add(NoteFromString(st, temps, longTemps,Aliases));
                            }
                        }
                        if (ligne != "")
                        {
                            temps += longTemps;
                        }
                    }
                }
                return retour;
            }
            private static Note NoteFromString(string NoteInfo, double temps, double longTemps, Dictionary<string,double> alias)
            {
                bool Anticiper = NoteInfo[0] == '@';
                if(Anticiper)
                {
                    NoteInfo = NoteInfo.Substring(1);
                }
                string[] splt1 = NoteInfo.Split('>');
                string infonote = splt1[0];
                bool Complete = infonote[infonote.Length - 1]=='@';
                if(Complete)
                {
                    infonote = infonote.Substring(0, infonote.Count() - 1);
                }
                string infoins = splt1[1];
                string[] splt2 = infonote.Split('/');
                double t = temps + calageFromString(splt2[2], longTemps, alias);
                int nt = HauteurFromString(splt2[0]);
                double duree = longueurFromString(splt2[1], longTemps, alias);
                double nuance;
                if(splt2.Count() >=4)
                {
                    nuance = nuanceFromString(splt2[3], alias);
                }
                else
                {
                    nuance = Math.Pow(10.0, (-40) / 30);
                }
                return new Note(t, duree, instrumentsFromString(infoins), nt,Anticiper,Complete,nuance);
            }
            private static int HauteurFromString(string note)
            {
                int octave = Convert.ToInt32(note.Split('.')[0]);
                int resultat = 12 * octave;
                string nt = note.Split('.')[1];
                switch (nt)
                {
                    case "do":
                        resultat += 0;
                        break;
                    case "do#":
                        resultat += 1;
                        break;
                    case "re":
                        resultat += 2;
                        break;
                    case "re#":
                        resultat += 3;
                        break;
                    case "mi":
                        resultat += 4;
                        break;
                    case "fa":
                        resultat += 5;
                        break;
                    case "fa#":
                        resultat += 6;
                        break;
                    case "sol":
                        resultat += 7;
                        break;
                    case "sol#":
                        resultat += 8;
                        break;
                    case "la":
                        resultat += 9;
                        break;
                    case "la#":
                        resultat += 10;
                        break;
                    case "si":
                        resultat += 11;
                        break;
                    default:
                        resultat += Convert.ToInt32(nt);
                        break;
                }
                return resultat;
            }
            private static double longueurFromString(string lg, double dureeTemps, Dictionary<string, double> alias)
            {
                double facteur = 0;
                switch(lg)
                {
                    case "cr":
                        facteur = 0.5;
                        break;
                    case "n":
                        facteur = 1;
                        break;
                    case "b":
                        facteur = 2;
                        break;
                    case "bp":
                        facteur = 4;
                        break;
                    default:
                        double delta;
                        if (double.TryParse(lg, out delta))
                        {
                            facteur = delta;
                        }
                        else
                        {
                            facteur += alias[lg];
                        }
                        break;
                }
                return dureeTemps * facteur;
            }
            private static double calageFromString(string cal,double dureeTemps, Dictionary<string, double> alias)
            {
                double facteur = 0;
                switch(cal)
                {
                    case "t":
                        facteur = 0;
                        break;
                    case "ct":
                        facteur = 0.5;
                        break;
                    default:
                        double delta;
                        if (double.TryParse(cal, out delta))
                        {
                            facteur = delta;
                        }
                        else
                        {
                            facteur = alias[cal];
                        }
                        break;
                }
                return facteur * dureeTemps;
            }
            private static double nuanceFromString(string nu, Dictionary<string, double> alias)
            {
                double facteur = 0;
                switch (nu)
                {
                    case "ppp":
                        facteur = 10;
                        break;
                    case "pp":
                        facteur = 20;
                        break;
                    case "p":
                        facteur = 30;
                        break;
                    case "mp":
                        facteur = 40;
                        break;
                    case "sotvo":
                        facteur = 10*(4+1.0/4.0);
                        break;
                    case "mezvo":
                        facteur = 10*(4+0.5);
                        break;
                    case "pocfor":
                        facteur = 10*(4+3.0/4.0);
                        break;
                    case "mf":
                        facteur = 50;
                        break;
                    case "f":
                        facteur = 60;
                        break;
                    case "ff":
                        facteur = 70;
                        break;
                    case "fff":
                        facteur = 80;
                        break;
                    default:
                        double delta;
                        if (double.TryParse(nu, out delta))
                        {
                            facteur = delta;
                        }
                        else
                        {
                            facteur = alias[nu];
                        }
                        break;
                }
                return Math.Pow(10.0,(facteur-80)/30);
            }
            private static List<int> instrumentsFromString(string ins)
            {
                List<int> retour = new List<int>();
                string[] inst = ins.Split('.');
                foreach(string instt in inst)
                {
                    retour.Add(Convert.ToInt32(instt));
                }
                return retour;
            }
        }
        class Orchestre
        {
            private Instrument[] Instruments;

            public Orchestre(Instrument[] instruments)
            {
                Instruments = instruments;
            }
            public Instrument GetInstrument(int n)
            {
                return Instruments[n];
            }
            public Instrument[] GetInstruments()
            {
                return Instruments;
            }
        }
        class Musique
        {
            private double fe;
            private Orchestre orchestre;
            private List<Note> notes;

            public Musique(double fe, Orchestre orchestre, List<Note> notes)
            {
                this.fe = fe;
                this.orchestre = orchestre;
                this.notes = notes;
            }

            public double[] Jouer()
            {
                List<double[]> sons = new List<double[]>();
                List<int> echantsDebuts = new List<int>();
                List<int> echantsFin = new List<int>();
                foreach (Note n in notes)
                {
                    double[] son = n.Jouer(orchestre, fe);
                    sons.Add(son);
                    int edeb = (int)(n.GetT0() * fe);
                    echantsDebuts.Add(edeb);
                    echantsFin.Add(edeb + son.GetLength(0));
                }
                int echMax = echantsFin.Max();
                double[] retour = new double[echMax];
                for(int i=0;i<echMax;i++)
                {
                    double somme = 0;
                    for(int j=0;j<sons.Count;j++)
                    {
                        if(i>=echantsDebuts.ElementAt(j) && i< echantsFin.ElementAt(j))
                        {
                            somme += sons.ElementAt(j)[i - echantsDebuts.ElementAt(j)];
                        }
                    }
                    retour[i] = somme;
                }
                double mx = retour.Max(i => Math.Abs(i));
                if(mx !=0)
                {
                    for(int i=0;i<echMax;i++)
                    {
                        retour[i] /= mx;
                    }
                }
                return retour;
            }
        }

    }
}
