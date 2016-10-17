#include "Quadrique.h"

using namespace Scene;

///////////////////////////////////////////////////////////////////////////////
///  public overloaded constructor  CQuadrique \n
///  Description : Constructeur par défaut
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CQuadrique::CQuadrique(void)
	: ISurface()
	, m_Quadratique(CVecteur3::ZERO)
	, m_Lineaire(CVecteur3::ZERO)
	, m_Mixte(CVecteur3::ZERO)
	, m_Cst(RENDRE_REEL(0))
{}

///////////////////////////////////////////////////////////////////////////////
///  public overloaded constructor  CQuadrique \n
///  Description : Constructeur par défaut
///
///  @param [in]       Quadric const Scene::CQuadrique &   la quadrique à copier
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CQuadrique::CQuadrique(const CQuadrique& Quadric)
	: ISurface(Quadric)
	, m_Quadratique(Quadric.m_Quadratique)
	, m_Lineaire(Quadric.m_Lineaire)
	, m_Mixte(Quadric.m_Mixte)
	, m_Cst(Quadric.m_Cst)
{}

///////////////////////////////////////////////////////////////////////////////
///  public virtual destructor  ~CQuadrique \n
///  Description : Destructeur
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CQuadrique::~CQuadrique(void)
{
}

///////////////////////////////////////////////////////////////////////////////
///  public  operator = \n
///  Description : Opérateur de copie
///
///  @param [in]       Quadric const Scene::CQuadrique &    La quadrique à copier
///
///  @return Scene::CQuadrique & La quadrique modifiée
///
///  @author Olivier Dionne 
///  @date   14/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CQuadrique& CQuadrique::operator = (const CQuadrique& Quadric)
{
	ISurface::operator =(Quadric);
	m_Quadratique = Quadric.m_Quadratique;
	m_Lineaire = Quadric.m_Lineaire;
	m_Mixte = Quadric.m_Mixte;
	m_Cst = Quadric.m_Cst;
	return (*this);
}

///////////////////////////////////////////////////////////////////////////////
///  protected virtual constant  AfficherInfoDebug \n
///  Description : Implémente le déboguage polymorphique par flux de sortie
///
///  @param [in, out]  Out std::ostream &    Le flux de sortie
///
///  @return std::ostream & Le flux de sortie modifié
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
ostream& CQuadrique::AfficherInfoDebug(ostream& Out) const
{
	Out << "[DEBUG]: Quadric.Quadratique       = " << m_Quadratique << endl;
	Out << "[DEBUG]: Quadric.Lineaire          = " << m_Lineaire << endl;
	Out << "[DEBUG]: Quadric.Mixte             = " << m_Mixte << endl;
	Out << "[DEBUG]: Quadric.Constante         = " << m_Cst;
	return Out;
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual  Pretraitement \n
///  Description : Pretraitement des données de la quadrique( appelé AVANT le lancer)
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
void CQuadrique::Pretraitement(void)
{
	// Algorithme tiré de ... 
	// R. Goldman, "Two Approach to a Computer Model for Quadric Surfaces",
	// IEEE CG&A, Sept 1983, pp.21

	REAL A = m_Quadratique.x;
	REAL B = m_Quadratique.y;
	REAL C = m_Quadratique.z;
	REAL D = m_Mixte.z    * RENDRE_REEL(0.5);
	REAL E = m_Mixte.x    * RENDRE_REEL(0.5);
	REAL F = m_Mixte.y    * RENDRE_REEL(0.5);
	REAL G = m_Lineaire.x * RENDRE_REEL(0.5);
	REAL H = m_Lineaire.y * RENDRE_REEL(0.5);
	REAL J = m_Lineaire.z * RENDRE_REEL(0.5);
	REAL K = m_Cst;

	CMatrice4 Q(A, D, F, G,
		D, B, E, H,
		F, E, C, J,
		G, H, J, K);

	CMatrice4 Inverse = m_Transformation.Inverse();

	Q = Inverse * Q * Inverse.Transpose();

	m_Quadratique.x = Q[0][0];
	m_Quadratique.y = Q[1][1];
	m_Quadratique.z = Q[2][2];
	m_Cst = Q[3][3];
	m_Mixte.x = Q[1][2] * RENDRE_REEL(2.0);
	m_Mixte.y = Q[0][2] * RENDRE_REEL(2.0);
	m_Mixte.z = Q[0][1] * RENDRE_REEL(2.0);
	m_Lineaire.x = Q[0][3] * RENDRE_REEL(2.0);
	m_Lineaire.y = Q[1][3] * RENDRE_REEL(2.0);
	m_Lineaire.z = Q[2][3] * RENDRE_REEL(2.0);
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual  Intersection \n
///  Description : Effectue l'intersection Rayon/Quadrique
///
///  @param [in]       Rayon const CRayon &    Le rayon à tester
///
///  @return Scene::CIntersection Le résultat de l'ntersection
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CIntersection CQuadrique::Intersection(const CRayon& Rayon)
{
	CIntersection Result;
	CVecteur3 dir = Rayon.ObtenirDirection();
	CVecteur3 orig = Rayon.ObtenirOrigine();

	//On calcule les coefficients de l'équation quadratique Aq.t²+Bq.t+Cq=0 obtenue en introduisant l'expression paramétrique de la droite dans l'équation de la quadrique
	REAL Aq = m_Quadratique.x * dir.x * dir.x
		+ m_Mixte.z * dir.x * dir.y
		+ m_Mixte.y * dir.x * dir.z
		+ m_Quadratique.y * dir.y * dir.y
		+ m_Mixte.x * dir.y * dir.z
		+ m_Quadratique.z * dir.z * dir.z;

	REAL Bq = 2 * m_Quadratique.x * dir.x * orig.x
		+ m_Mixte.z*(orig.x*dir.y + orig.y*dir.x)
		+ m_Mixte.y*(orig.x*dir.z + orig.z*dir.x)
		+ m_Lineaire.x*dir.x
		+ 2 * m_Quadratique.y*dir.y*orig.y
		+ m_Mixte.x*(dir.z*orig.y + dir.y*orig.z)
		+ m_Lineaire.y*dir.y
		+ 2 * m_Quadratique.z*dir.z*orig.z
		+ m_Lineaire.z*dir.z;

	REAL Cq = m_Quadratique.x*orig.x*orig.x
		+ m_Mixte.z*orig.x*orig.y
		+ m_Mixte.y*orig.x*orig.z
		+ m_Lineaire.x*orig.x
		+ m_Quadratique.y*orig.y*orig.y
		+ m_Mixte.x*orig.y*orig.z
		+ m_Lineaire.y*orig.y
		+ m_Quadratique.z*orig.z*orig.z
		+ m_Lineaire.z*orig.z
		+ m_Cst;

	// valeur par défaut: t=0 correspond à l'origine de la droite, pas d'intersection possible
	REAL t = 0;
	if (Aq == 0 && Bq != 0)
		//Une solution unique
		t = -Cq / Bq;
	else {
		//Déterminant
		REAL delta = Bq*Bq - 4 * Aq*Cq;
		if (delta >= 0) {
			//Deux solutions (éventuellement égales)
			REAL t0 = (-Bq - sqrt(delta)) / (2 * Aq);
			REAL t1 = (-Bq + sqrt(delta)) / (2 * Aq);
			//On choisit le point le plus proche de la caméra
			t = t0 < t1 ? t0 : t1;
		}
	}

	// En cas d'intersection
	if (t > 0) {
		Result.AjusterSurface(this);
		Result.AjusterDistance(CVecteur3::Norme(dir*t));

		//Coordonnées du point d'intersection
		REAL x = orig.x + dir.x * t;
		REAL y = orig.y + dir.y * t;
		REAL z = orig.z + dir.z * t;

		//Normale à la quadrique au point d'intersection
		REAL xn = 2 * m_Quadratique.x*x
			+ m_Mixte.z * y
			+ m_Mixte.y * z
			+ m_Lineaire.x;
		REAL yn = m_Mixte.z * x
			+ 2 * m_Quadratique.y*y
			+ m_Mixte.x*z
			+ m_Lineaire.y;
		REAL zn = m_Mixte.y *x
			+ m_Mixte.x*y
			+ 2 * m_Quadratique.z*z
			+ m_Lineaire.z;

		Result.AjusterNormale(CVecteur3::Normaliser(CVecteur3(xn, yn, zn)));
	}

	return Result;
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual constant  Copier \n
///  Description : Alloue une copie de la quadrique courante
///
///  @return Scene::CQuadrique *la copie de la quadrique
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CQuadrique* CQuadrique::Copier(void) const
{
	return new CQuadrique(*this);
}