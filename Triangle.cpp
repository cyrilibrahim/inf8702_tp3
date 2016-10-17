#include "Triangle.h"

using namespace Scene;

///////////////////////////////////////////////////////////////////////////////
///  public overloaded constructor  CTriangle \n
///  Description : Constructeur par défaut
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CTriangle::CTriangle( void )
    : ISurface (                 )
    , m_Normale( CVecteur3::ZERO )
{}

///////////////////////////////////////////////////////////////////////////////
///  public overloaded constructor  CTriangle \n
///  Description : Constructeur par défaut
///
///  @param [in]       Triangle const Scene::CTriangle &    Le triangle à copier
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CTriangle::CTriangle( const CTriangle& Triangle )
    : ISurface ( Triangle           )
    , m_Normale( Triangle.m_Normale )
{
    for( int i = 0; i < 3; i++ )
        m_Pts[ i ] = Triangle.m_Pts[ i ];
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual destructor  ~CTriangle \n
///  Description : Destructeur
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CTriangle::~CTriangle( void )
{
}

///////////////////////////////////////////////////////////////////////////////
///  public  operator = \n
///  Description : Opérateur de copie
///
///  @param [in]       Triangle const Scene::CTriangle &    Le triangle à copier
///
///  @return Scene::CTriangle & Le triangle modifié
///
///  @author Olivier Dionne 
///  @date   14/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CTriangle& CTriangle::operator = ( const CTriangle& Triangle )
{
    ISurface::operator =( Triangle );
    m_Normale = Triangle.m_Normale;

    for( int i = 0; i < 3; i++ )
        m_Pts[ i ] = Triangle.m_Pts[ i ];

    return ( *this );
}

///////////////////////////////////////////////////////////////////////////////
///  protected virtual constant  AfficherInfoDebug \n
///  Description : Implémente le déboguage polymorphique par flux de sortie
///
///  @param [in, out]  Out std::ostream &   Le flux de sortie
///
///  @return std::ostream & Le flux de sortie modifié
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
ostream& CTriangle::AfficherInfoDebug( ostream& Out ) const
{
    Out << "[DEBUG]: Triangle.Point1         = " << m_Pts[ 0 ] << endl;
    Out << "[DEBUG]: Triangle.Point2         = " << m_Pts[ 1 ] << endl;
    Out << "[DEBUG]: Triangle.Point3         = " << m_Pts[ 2 ] << endl;
    Out << "[DEBUG]: Triangle.Normale        = " << m_Normale;
    return Out;
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual  Pretraitement \n
///  Description : Pretraitement des données du triangle ( Appelé AVANT le lancer)
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
void CTriangle::Pretraitement( void )
{
    for( int i = 0; i < 3; i++ )
        m_Pts[ i ] = m_Pts[ i ] * m_Transformation;
    CalculerNormale();
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual  Intersection \n
///  Description : Effectue l'intersection Rayon/Triangle
///
///  @param [in]       Rayon const CRayon &    Le rayon à tester
///
///  @return Scene::CIntersection Le résultat de l'ntersection
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CIntersection CTriangle::Intersection( const CRayon& Rayon )
{
	CIntersection Result;

	// À COMPLÉTER ... 

	// Voici deux références pour acomplir le développement :
	// 1) Tomas Akenine-Moller and Eric Haines "Real-Time Rendering 2nd Ed." 2002, p.581
	// 2) Son article: http://www.graphics.cornell.edu/pubs/1997/MT97.pdf

	// Notez que la normale du triangle est déjà calculée lors du prétraitement
	// il suffit que de la passer à la structure d'intersection.
	
	CVecteur3 dirRayon = Rayon.ObtenirDirection();
	CVecteur3 origineRayon = Rayon.ObtenirOrigine();

	//E1 dans l'équation (edge 1)
	CVecteur3 Edge1 = m_Pts[1] - m_Pts[0];
	//E2 dans l'équation (edge 2)
	CVecteur3 Edge2 = m_Pts[2] - m_Pts[0];

	//Cacule de la variable P (pour le calcul du determinant) de la formule
	CVecteur3 P = CVecteur3::ProdVect(dirRayon, Edge2);

	//Calcul du determinant 
	REAL determinant = CVecteur3::ProdScal(P, Edge1);

	//Si le determinant est nul alors il n'y a pas d'intersection
	//On fait la comparaison à EPSILON au lieu de 0 pour des questions de precision
	if (determinant > -EPSILON && determinant < EPSILON)
		return Result;

	//Calcul de la variable T  (distance du sommet 0 a l'origine du rayon)
	CVecteur3 T = origineRayon - m_Pts[0];

	//Calcul de Q (variable intermediaire de calcul)
	CVecteur3 Q = CVecteur3::ProdVect(T, Edge1);

	//On calcul les valeur de u et v (coordonnées barycentrique de l'intersection sur le triangle)
	REAL u = (1 / determinant) * CVecteur3::ProdScal(P, T);
	
	//Si u est inferieur à 0 ou superieur à 1 cela signifie que l'intersection n'est pas dans le triangle
	if (u < 0 || u > 1)
		return Result;
	
	REAL v = (1 / determinant) * CVecteur3::ProdScal(Q, dirRayon);

	//Si v est inferieur à 0 ou v + u superieur à 1 cela signifie que l'intersection n'est pas dans le triangle
	if (v < 0 || (v + u) > 1)
		return Result;

	//Calcul du t à l'intersection de la formule du rayon (R(t) = Ro + Rd*t) 
	REAL t = (1 / determinant) * CVecteur3::ProdScal(Q, Edge2);

	//Calcul de la distance entre l'origine du rayon et le point d'interesction
	REAL distance = CVecteur3::Norme( t * dirRayon );

	//On ajuste le resultat d'intersection
	Result.AjusterSurface(this);
	Result.AjusterDistance(distance);
	Result.AjusterNormale(m_Normale);

    return Result;
}

///////////////////////////////////////////////////////////////////////////////
///  public virtual constant  Copier \n
///  Description : Alloue une copie du triangle courant
///
///  @return Scene::CTriangle * Nouvelle copie du triangle 
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CTriangle* CTriangle::Copier( void ) const
{
    return new CTriangle( *this );
}

///////////////////////////////////////////////////////////////////////////////
///  private  CalculerNormale \n
///  Description : Calculer la normale du triangle à partir des côtés
///
///  @return None
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
void CTriangle::CalculerNormale( void )
{
    CVecteur3 Edge1 = m_Pts[ 1 ] - m_Pts[ 0 ];
    CVecteur3 Edge2 = m_Pts[ 2 ] - m_Pts[ 0 ];
    m_Normale = CVecteur3::Normaliser( CVecteur3::ProdVect( Edge1, Edge2 ) );
}
