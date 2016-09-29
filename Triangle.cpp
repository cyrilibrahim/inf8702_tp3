#include "Triangle.h"

using namespace Scene;

///////////////////////////////////////////////////////////////////////////////
///  public overloaded constructor  CTriangle \n
///  Description : Constructeur par d�faut
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
///  Description : Constructeur par d�faut
///
///  @param [in]       Triangle const Scene::CTriangle &    Le triangle � copier
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
///  Description : Op�rateur de copie
///
///  @param [in]       Triangle const Scene::CTriangle &    Le triangle � copier
///
///  @return Scene::CTriangle & Le triangle modifi�
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
///  Description : Impl�mente le d�boguage polymorphique par flux de sortie
///
///  @param [in, out]  Out std::ostream &   Le flux de sortie
///
///  @return std::ostream & Le flux de sortie modifi�
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
///  Description : Pretraitement des donn�es du triangle ( Appel� AVANT le lancer)
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
///  @param [in]       Rayon const CRayon &    Le rayon � tester
///
///  @return Scene::CIntersection Le r�sultat de l'ntersection
///
///  @author Olivier Dionne 
///  @date   13/08/2008
///
///////////////////////////////////////////////////////////////////////////////
CIntersection CTriangle::Intersection( const CRayon& Rayon )
{
	CIntersection Result;

	//Ar�tes du triangle
	CVecteur3 e1 = m_Pts[1]-m_Pts[0];
	CVecteur3 e2 = m_Pts[2]-m_Pts[0];
	//Expressions r�utilis�es pour calculer la solution
	CVecteur3 T = Rayon.ObtenirOrigine() - m_Pts[0];
	CVecteur3 P = CVecteur3::ProdVect(Rayon.ObtenirDirection(),e2);
	CVecteur3 Q = CVecteur3::ProdVect(T,e1);

	//[t,u,v] : param�tres au point d'intersection
	CVecteur3 tuv = (1/CVecteur3::ProdScal(P,e1))*CVecteur3(
			CVecteur3::ProdScal(Q,e2),
			CVecteur3::ProdScal(P,T),
			CVecteur3::ProdScal(Q,Rayon.ObtenirDirection()));

	//Le point est � l'int�rieur du triangle si u>=0, v>=0, u+v<=1
	if (tuv.y>=0 && tuv.z>=0 && tuv.y+tuv.z<=1) {
		Result.AjusterSurface(this);
		Result.AjusterDistance(CVecteur3::Norme(Rayon.ObtenirDirection()*tuv.x));
		Result.AjusterNormale(m_Normale);
	}

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
///  Description : Calculer la normale du triangle � partir des c�t�s
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
