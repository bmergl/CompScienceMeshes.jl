
abstract type HexahedralElement{P} end


# -----------------------------------------------------------------------------
# Trilinearly mapped Hexahedron (standard)
# -----------------------------------------------------------------------------

struct Hexahedron{P} <: HexahedralElement{P}
    vertices::SVector{8, P}
    p21::P   # p2 - p1
    p41::P   # p4 - p1
    p51::P   # p5 - p1
    p1243::P # p1 - p2 - p4 + p3
    p1256::P # p1 - p2 - p5 + p6
    p1458::P # p1 - p4 - p5 + p8 
    p1to8::P #-p1 + p2 - p3 + p4 + p5 - p6 + p7 - p8
end
Hexahedron(p1::P, p2::P, p3::P, p4::P, p5::P, p6::P, p7::P, p8::P) where P = Hexahedron{P}(
    SVector(p1, p2, p3, p4, p5, p6, p7, p8), 
    p2 - p1,
    p4 - p1, 
    p5 - p1,
    p1 - p2 - p4 + p3,
    p1 - p2 - p5 + p6,
    p1 - p4 - p5 + p8,
    -p1 + p2 - p3 + p4 + p5 - p6 + p7 - p8
)
Hexahedron(p) = Hexahedron(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])

function coordtype(hex::Hexahedron{P}) where P eltype(P) end

function vertices(hex::Hexahedron) hex.vertices end

function cartesian(hex::Hexahedron, u)
    u1, u2, u3 = u

    return hex.vertices[1] + u1*hex.p21 + u2*hex.p41 + u3*hex.p51 + 
    u1*u2*hex.p1243 +
    u1*u3*hex.p1256 +
    u2*u3*hex.p1458 +
    u1*u2*u3*hex.p1to8
end

function tangents(hex::Hexahedron, u)
    u1, u2, u3 = u

    ∂ru1 = hex.p21 + u2*hex.p1243 + u3*hex.p1256 + u2*u3*hex.p1to8
    ∂ru2 = hex.p41 + u1*hex.p1243 + u3*hex.p1458 + u1*u3*hex.p1to8
    ∂ru3 = hex.p51 + u1*hex.p1256 + u2*hex.p1458 + u1*u2*hex.p1to8

    #return ∂ru1, ∂ru2, ∂ru3

    return SMatrix{3,3}(
        ∂ru1[1], ∂ru1[2], ∂ru1[3], # ---> column 1 for an SMatrix
        ∂ru2[1], ∂ru2[2], ∂ru2[3], # ---> column 2 for an SMatrix
        ∂ru3[1], ∂ru3[2], ∂ru3[3]  # ---> column 3 for an SMatrix
    )
end

function jacobian(hex::Hexahedron, u)
    J = tangents(hex, u)
    return dot(cross(J[:,1], J[:,2]), J[:,3])
end

function jacobian_(hex::Hexahedron, ∂ru1, ∂ru2, ∂ru3)
    return dot(cross(∂ru1, ∂ru2), ∂ru3)
end


# -----------------------------------------------------------------------------
# Affinely mapped hexahedron 
# -----------------------------------------------------------------------------

struct AffineHexahedron{P,T}  <: HexahedralElement{P}
    p1::P
    a::P   # p2 - p1
    b::P   # p4 - p1
    c::P   # p5 - p1
    volume::T
end
function AffineHexahedron(p1::P, a::P, b::P, c::P) where P 
    return AffineHexahedron{P,eltype(P)}(p1, a, b, c, dot(cross(a, b), c))
end
function AffineHexahedron(p1::P, p2::P, p3::P, p4::P, p5::P, p6::P, p7::P, p8::P) where P
    a = p2 - p1 # No test to see if it is really affinely mapped!
    b = p4 - p1 
    c = p5 - p1
    return AffineHexahedron(p1, a, b, c)
end

function coordtype(hex::AffineHexahedron{P}) where P eltype(P) end

function vertices(hex::AffineHexahedron{P}) where P 
    a = hex.a
    b = hex.b
    c = hex.c

    p1 = hex.p1
    p2 = p1 + a
    p3 = p2 + b
    p4 = p1 + b

    p5 = p1 + c
    p6 = p5 + a
    p7 = p6 + b
    p8 = p5 + b

    return SVector(p1, p2, p3, p4, p5, p6, p7, p8)
end

function cartesian(hex::AffineHexahedron, u)
    return hex.p1 + u[1]*hex.a + u[2]*hex.b + u[3]*hex.c
end

function tangents(hex::AffineHexahedron, u)

    ∂ru1 = hex.a
    ∂ru2 = hex.b
    ∂ru3 = hex.c

    return SMatrix{3,3}(
        ∂ru1[1], ∂ru1[2], ∂ru1[3], # ---> column 1 for an SMatrix
        ∂ru2[1], ∂ru2[2], ∂ru2[3], # ---> column 2 for an SMatrix
        ∂ru3[1], ∂ru3[2], ∂ru3[3]  # ---> column 3 for an SMatrix
    )
end

function jacobian(hex::AffineHexahedron, u)
    return hex.volume
end


# -----------------------------------------------------------------------------
#  RefHexahedron (for dispatch), RefHexahedron_ (for fast computation)
# -----------------------------------------------------------------------------

struct RefHexahedron_{P} <: HexahedralElement{P}
    p1::P
    a::P # p2-p1
    b::P # p4-p1
    c::P # p5-p1
end

function RefHexahedron_(p1::P, p2::P, p3::P, p4::P, p5::P, p6::P, p7::P, p8::P) where P
    return RefHexahedron_(p1, p2-p1, p4-p1, p5-p1)
end
function cartesian(hex::RefHexahedron_, u)
    return hex.p1 + u[1]*hex.a + u[2]*hex.b + u[3]*hex.c
end


struct RefHexahedron{T} end

domain(hex::HexahedralElement{P}) where {P} = RefHexahedron{eltype(P)}()

function vertices(hex::RefHexahedron{T}) where T
    SVector(
        point(T, 0, 0, 0), #p1 , A
        point(T, 1, 0, 0), #p2 , B
        point(T, 1, 1, 0), #p3 , C
        point(T, 0, 1, 0), #p4 , D
        point(T, 0, 0, 1), #p5 , E
        point(T, 1, 0, 1), #p6 , F
        point(T, 1, 1, 1), #p7 , G
        point(T, 0, 1, 1)  #p8 , H
        )
end

function permute_vertices(hex::RefHexahedron, I)
    V = vertices(hex)
    return AffineHexahedron(V[I[1]], V[I[2]], V[I[3]], V[I[4]], V[I[5]], V[I[6]], V[I[7]], V[I[8]])
end


# -----------------------------------------------------------------------------
# Neighborhood
# -----------------------------------------------------------------------------

struct NeighborhoodHex{C,P,Q,T,J}
    chart::C
    parametric::P
    cartesian::Q
    tangents::T
    jacobian::J
end

function neighborhood(hex::Hexahedron, u)
    c = cartesian(hex, u)
    J = tangents(hex, u)
    j = jacobian_(hex, J[:,1], J[:,2], J[:,3])
    return NeighborhoodHex(hex, u, c, J, j)
end
function neighborhood(hex::AffineHexahedron, u)
    c = cartesian(hex, u)
    J = tangents(hex, u)
    return NeighborhoodHex(hex, u, c, J, hex.volume)
end

function parametric(mp::NeighborhoodHex) mp.parametric end
function cartesian(mp::NeighborhoodHex) mp.cartesian end
function tangents(mp::NeighborhoodHex) mp.tangents end
function tangents(mp::NeighborhoodHex, i::Int) mp.tangents[:,i] end
function jacobian(mp::NeighborhoodHex) mp.jacobian end


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

@testitem "hexahedron" begin

    # cube
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    p5 = point(0,0,1)
    p6 = point(1,0,1)
    p7 = point(1,1,1)
    p8 = point(0,1,1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)

    @test coordtype(hex) == Float64

    mp = neighborhood(hex, (0.5,0.5,0.5))
    @test cartesian(mp) ≈ point(0.5,0.5,0.5)
    @test tangents(mp, 1) ≈ point(1,0,0)
    @test tangents(mp, 2) ≈ point(0,1,0)
    @test tangents(mp, 3) ≈ point(0,0,1)
    @test jacobian(hex, (0.5,0.5,0.5)) ≈ 1.0
    @test hex.p1243 ≈ point(0,0,0)
    @test hex.p1256 ≈ point(0,0,0)
    @test hex.p1458 ≈ point(0,0,0)
    @test hex.p1to8 ≈ point(0,0,0)


    # distorted hexahedron
    p1 = point(0.0, 0.0, 0.0)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)
    p5 = point(0.0, 0.0, 0.9)
    p6 = point(2.3, -0.2, 5.2)
    p7 = point(2.7, 1.1, 3.6)
    p8 = point(-0.2, 2.2, 3.1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)

    @test cartesian(neighborhood(hex, (0,0,0))) ≈ p1
    @test cartesian(neighborhood(hex, (1,0,0))) ≈ p2
    @test cartesian(neighborhood(hex, (1,1,0))) ≈ p3 
    @test cartesian(neighborhood(hex, (0,1,0))) ≈ p4 
    @test cartesian(neighborhood(hex, (0,0,1))) ≈ p5 
    @test cartesian(neighborhood(hex, (1,0,1))) ≈ p6
    @test cartesian(neighborhood(hex, (1,1,1))) ≈ p7 
    @test cartesian(neighborhood(hex, (0,1,1))) ≈ p8
    mp = neighborhood(hex, (0,0,0))
    @test tangents(mp, 1) ≈ hex.p21
    @test tangents(mp, 2) ≈ hex.p41
    @test tangents(mp, 3) ≈ hex.p51
    @test jacobian(hex, (0.5,0.5,0.5)) > 0.0
    @test cartesian(neighborhood(hex, (0.5,0.5,0.5))) ≈ sum(hex.vertices)/8
end



@testitem "affine hexahedron" begin
    # cube
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    p5 = point(0,0,1)
    p6 = point(1,0,1)
    p7 = point(1,1,1)
    p8 = point(0,1,1)

    hex = CompScienceMeshes.AffineHexahedron(p1,p2,p3,p4,p5,p6,p7,p8)

    @test coordtype(hex) == Float64

    mp = neighborhood(hex, (0.5,0.5,0.5))
    @test cartesian(mp) ≈ point(0.5,0.5,0.5)
    @test tangents(mp, 1) ≈ point(1,0,0)
    @test tangents(mp, 2) ≈ point(0,1,0)
    @test tangents(mp, 3) ≈ point(0,0,1)
    @test hex.volume ≈ jacobian(hex, (0.5,0.5,0.5))
    @test hex.volume ≈ 1.0
end


@testitem "reference hexahedron" begin

    # distorted hexahedron
    p1 = point(0.0, 0.0, 0.0)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)
    p5 = point(0.0, 0.0, 0.9)
    p6 = point(2.3, -0.2, 5.2)
    p7 = point(2.7, 1.1, 3.6)
    p8 = point(-0.2, 2.2, 3.1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)
    
    refhex = domain(hex)
    @test typeof(refhex).parameters[1] == Float64
    @test vertices(refhex)[5] == point(0, 0, 1)


    I = [4,1,2,3,8,5,6,7]
    refhex_ = CompScienceMeshes.permute_vertices(refhex, I)

    @test refhex_.p1 ≈ point(0.0, 1.0, 0.0)
    @test refhex_.a ≈ point(0.0, 0.0, 0.0) - point(0.0, 1.0, 0.0)
    @test refhex_.b ≈ point(1.0, 1.0, 0.0) - point(0.0, 1.0, 0.0)
    @test refhex_.c ≈ point(0.0, 0.0, 1.0) 

    @test cartesian(refhex_, (0.3,0.1,0.6)) ≈ refhex_.p1 + 0.3*refhex_.a + 0.1*refhex_.b + 0.6*refhex_.c
end












# -----------------------------------------------------------------------------
# QuadrilateralElement from HexahedralElement
# -----------------------------------------------------------------------------

struct QuadFromHex{P,T,Q,H}  <: CompScienceMeshes.QuadrilateralElement{P}
    index::Int # 1 to 6
    quad::Q  #QuadrilateralElement
    hex::H # HexahedralElement
end

function coordtype(q::QuadFromHex{P,T}) where {P,T} T end
function vertices(q::QuadFromHex) vertices(q.quad) end
function cartesian(q::QuadFromHex, u) cartesian(q.quad, u) end
function tangents(q::QuadFromHex, u) tangents(q.quad,u) end
function normal(q::QuadFromHex, u) normal(q.quad, u) end
function jacobian(q::QuadFromHex, u) jacobian(q.quad) end

function cellboundaryfacets(hex::Hexahedron{P}) where P
    p = vertices(hex)

    quad1 = Quadrilateral(p[1],p[4],p[3],p[2])
    quad2 = Quadrilateral(p[5],p[6],p[7],p[8])
    quad3 = Quadrilateral(p[1],p[2],p[6],p[5])
    quad4 = Quadrilateral(p[3],p[4],p[8],p[7])
    quad5 = Quadrilateral(p[1],p[5],p[8],p[4])
    quad6 = Quadrilateral(p[2],p[3],p[7],p[6])

    T = coordtype(hex)
    Q = typeof(quad1)
    H = typeof(hex)

    return (
        QuadFromHex{P,T,Q,H}(1,quad1,hex), 
        QuadFromHex{P,T,Q,H}(2,quad2,hex), 
        QuadFromHex{P,T,Q,H}(3,quad3,hex), 
        QuadFromHex{P,T,Q,H}(4,quad4,hex), 
        QuadFromHex{P,T,Q,H}(5,quad5,hex), 
        QuadFromHex{P,T,Q,H}(6,quad6,hex)
    )
end

function cellboundaryfacets(hex::AffineHexahedron{P}) where P
    p = vertices(hex)

    quad1 = AffineQuadrilateral(p[1],p[4],p[3],p[2])
    quad2 = AffineQuadrilateral(p[5],p[6],p[7],p[8])
    quad3 = AffineQuadrilateral(p[1],p[2],p[6],p[5])
    quad4 = AffineQuadrilateral(p[3],p[4],p[8],p[7])
    quad5 = AffineQuadrilateral(p[1],p[5],p[8],p[4])
    quad6 = AffineQuadrilateral(p[2],p[3],p[7],p[6])

    T = coordtype(hex)
    Q = typeof(quad1)
    H = typeof(hex)

    return (
        QuadFromHex{P,T,Q,H}(1,quad1,hex), 
        QuadFromHex{P,T,Q,H}(2,quad2,hex), 
        QuadFromHex{P,T,Q,H}(3,quad3,hex), 
        QuadFromHex{P,T,Q,H}(4,quad4,hex), 
        QuadFromHex{P,T,Q,H}(5,quad5,hex), 
        QuadFromHex{P,T,Q,H}(6,quad6,hex)
    )
end

function _quad2hexcoords(quadfromhex::QuadFromHex{P,T}, u) where {P,T}
    index = quadfromhex.index
    index == 1 && return (u[2], u[1], T(0))  # Quadrilateral(p1,p2,p3,p4)
    index == 2 && return (u[1], u[2], T(1))  # Quadrilateral(p5,p6,p7,p8)
    index == 3 && return (u[1], T(0),  u[2]) # Quadrilateral(p1,p2,p6,p5) 
    index == 4 && return (T(1)-u[1], T(1),  u[2]) # Quadrilateral(p3,p4,p8,p7)
    index == 5 && return (T(0), u[2], u[1])  # Quadrilateral(p1,p5,p8,p4)
    index == 6 && return (T(1), u[1], u[2])  # Quadrilateral(p2,p3,p7,p6)

    error("index=$(index), index=1...6 is allowed.")
end


# -----------------------------------------------------------------------------
# Neighborhood: QuadrilateralElement from HexahedralElement
# -----------------------------------------------------------------------------

struct NeighborhoodQuadFromHex{N1,N2}
    nbquad::N1 
    nbhex::N2
end

function neighborhood(el::QuadFromHex, u) # u=(u1,u2) vom Quadrilateral
    nb_quad = neighborhood(el.quad, u)
    u_hex = _quad2hexcoords(el, u)
    nb_hex = neighborhood(el.hex, u_hex)
    return NeighborhoodQuadFromHex(nb_quad,nb_hex)
end

function parametric(nb::NeighborhoodQuadFromHex) nb.nbquad.parametric end
function cartesian(nb::NeighborhoodQuadFromHex) nb.nbquad.cartesian end
function tangents(nb::NeighborhoodQuadFromHex) nb.nbquad.tangents end
function tangents(nb::NeighborhoodQuadFromHex, i::Int) nb.nbquad.tangents[:,i] end
function jacobian(nb::NeighborhoodQuadFromHex) nb.nbquad.jacobian end
function normal(nb::NeighborhoodQuadFromHex) nb.nbquad.normal end
function parent_neighboorhood(nb::NeighborhoodQuadFromHex) nb.nbhex end


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

@testitem "QuadFromHex" begin

    using LinearAlgebra

    # distorted hexahedron
    p1 = point(0.0, 0.0, 0.0)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)
    p5 = point(0.0, 0.0, 0.9)
    p6 = point(2.3, -0.2, 5.2)
    p7 = point(2.7, 1.1, 3.6)
    p8 = point(-0.2, 2.2, 3.1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)
    quads = CompScienceMeshes.cellboundaryfacets(hex)

    # test _quad2hexcoords
    for i in 1:6
        nb = neighborhood(quads[i], (0.118,0.609))
        nb_ = CompScienceMeshes.parent_neighboorhood(nb)
        @test cartesian(nb) ≈ cartesian(nb_)
    end

    # test normal
    for i in 1:6
        nb = neighborhood(quads[i], (0.5,0.5))
        center = cartesian(nb)
        n = normal(nb)

        nb0 = neighborhood(hex, (0.5,0.5,0.5))
        center0 = cartesian(nb0)

        outward = center - center0

        @test dot(n,outward) > 0.0
    end


    # cube
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    p5 = point(0,0,1)
    p6 = point(1,0,1)
    p7 = point(1,1,1)
    p8 = point(0,1,1)

    hex = CompScienceMeshes.AffineHexahedron(p1,p2,p3,p4,p5,p6,p7,p8)
    quads = CompScienceMeshes.cellboundaryfacets(hex)

    # test _quad2hexcoords
    for i in 1:6
        nb = neighborhood(quads[i], (0.118,0.609))
        nb_ = CompScienceMeshes.parent_neighboorhood(nb)
        @test cartesian(nb) ≈ cartesian(nb_)
    end

    # test normal
    for i in 1:6
        nb = neighborhood(quads[i], (0.5,0.5))
        center = cartesian(nb)
        n = normal(nb)

        nb0 = neighborhood(hex, (0.5,0.5,0.5))
        center0 = cartesian(nb0)

        outward = center - center0

        @test dot(n, outward) > 0.0
    end

end


