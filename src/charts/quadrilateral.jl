
abstract type QuadrilateralElement{P} end


# -----------------------------------------------------------------------------
# Bilinearly mapped Quadrilateral 
# -----------------------------------------------------------------------------

struct Quadrilateral{P} <: QuadrilateralElement{P}
    vertices::SVector{4, P}
    p21::P   # p2 - p1
    p41::P   # p4 - p1
    p1243::P # p1 - p2 - p4 + p3 
end
function Quadrilateral(p1::P, p2::P, p3::P, p4::P) where P
    return Quadrilateral{P}(SVector(p1,p2,p3,p4), p2-p1, p4-p1, p1-p2-p4+p3)
end
function Quadrilateral(p)
    @assert length(p) == 4
    return Quadrilateral(p[1], p[2], p[3], p[4])
end


function coordtype(quad::Quadrilateral{P}) where {P} eltype(P) end

function vertices(quad::Quadrilateral) quad.vertices end

function cartesian(quad::Quadrilateral, u)
    u1, u2 = u
    return quad.vertices[1] + u1*quad.p21 + u2*quad.p41 + u1*u2*quad.p1243
end

function tangents(quad::Quadrilateral, u)
    ∂ru = quad.p21 + u[2] * quad.p1243
    ∂rv = quad.p41 + u[1] * quad.p1243
    return hcat(∂ru, ∂rv)
end

function normal(quad::Quadrilateral, u)
    Du = tangents(quad, u)
    return normalize(Du[:,1] × Du[:,2])
end

function jacobian(quad::Quadrilateral, u)
    Du = tangents(quad, u)
    g = Du' * Du
    return √det(g)
end

function faces(ch::Quadrilateral)
    p = vertices(ch)
    return SVector(
        simplex(p[1],p[2]),
        simplex(p[2],p[3]),
        simplex(p[3],p[4]),
        simplex(p[4],p[1]))
end


# -----------------------------------------------------------------------------
#  Affinely mapped Quadrilateral 
# -----------------------------------------------------------------------------

struct AffineQuadrilateral{P,T} <: QuadrilateralElement{P}
    p1::P
    a::P   # p2 - p1
    b::P   # p4 - p1
    n::P # normal
    volume::T
end

function AffineQuadrilateral(p1::P, a::P, b::P) where P 
    cross_a_b = cross(a, b)  

    return AffineQuadrilateral{P,eltype(P)}(p1, a, b, normalize(cross_a_b), norm(cross_a_b))
end
function AffineQuadrilateral(p1::P, p2::P, p3::P, p4::P) where P
    return AffineQuadrilateral(p1, p2 - p1, p4 - p1) # No test to see if it is really affinely mapped!
end


function coordtype(quad::AffineQuadrilateral{P}) where P eltype(P) end

function vertices(quad::AffineQuadrilateral{P}) where P 
    a = quad.a
    b = quad.b

    p1 = quad.p1
    p2 = p1 + a
    p3 = p2 + b
    p4 = p1 + b

    return SVector(p1, p2, p3, p4)
end

function cartesian(quad::AffineQuadrilateral, u)
    return quad.p1 + u[1]*quad.a + u[2]*quad.b
end

function tangents(quad::AffineQuadrilateral, u)
    ∂ru = quad.a
    ∂rv = quad.b
    return hcat(∂ru, ∂rv)
end

function jacobian(quad::AffineQuadrilateral, u)
    return quad.volume
end


# -----------------------------------------------------------------------------
#  RefQuadrilateral (for dispatch), RefQuadrilateral_ (for fast computation)
# -----------------------------------------------------------------------------

struct RefQuadrilateral_{P} <: QuadrilateralElement{P}
    p1::P # P=SVector{2,T} intended!
    a::P # p2-p1
    b::P # p4-p1
end
function RefQuadrilateral_(p1::P, p2::P, p3::P, p4::P) where P
    return RefQuadrilateral_{P}(p1, p2-p1, p4-p1)
end
function cartesian(quad::RefQuadrilateral_, u) #!!!!
    return quad.p1 + u[1]*quad.a + u[2]*quad.b # 2D points
end


struct RefQuadrilateral{T} end

function domain(quad::QuadrilateralElement{P}) where P return RefQuadrilateral{eltype(P)}() end

function vertices(ch::RefQuadrilateral{T}) where {T}
    SVector(
        point(T,0,0),
        point(T,1,0),
        point(T,1,1),
        point(T,0,1))
end

neighborhood(quad::RefQuadrilateral, u) = SVector(u)

neighborhood(ch::RefQuadrilateral, u::AbstractVector) = SVector{length(u)}(u)

function permute_vertices(ch::RefQuadrilateral, I)
    V = vertices(ch)
    return RefQuadrilateral_(V[I[1]], V[I[2]], V[I[3]], V[I[4]])
end

function faces(ch::RefQuadrilateral)
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    return SVector(
        simplex(p1,p2),
        simplex(p2,p3),
        simplex(p3,p4),
        simplex(p4,p1),)
end

function quadpoints(ch::RefQuadrilateral{T}, rule) where {T}

    U1, W1 = legendre(rule, zero(T), one(T))
    U2, W2 = legendre(rule, zero(T), one(T))

    [(neighborhood(ch, (u1,u2)), w1*w2) for (u1,w1) in zip(U1,W1) for (u2,w2) in zip(U2,W2)]
end

function permute_vertices(q::Quadrilateral, I)
    verts = vertices(q)[I]
    return Quadrilateral(verts[1], verts[2], verts[3], verts[4])
end

function permute_vertices(q::AffineQuadrilateral, I)
    verts = vertices(q)[I]
    return AffineQuadrilateral(verts[1], verts[2], verts[3], verts[4])
end

function center(q::QuadrilateralElement)
    T = coordtype(q)
    h = T(0.5)
    return neighborhood(q, (h,h))
end



# -----------------------------------------------------------------------------
#  Neighborhood
# -----------------------------------------------------------------------------

struct NeighborhoodQuad{C,P,Q,T,J,N}
    chart::C
    parametric::P
    cartesian::Q
    tangents::T
    jacobian::J
    normal::N
end
function neighborhood(quad::Quadrilateral, u)
    c = cartesian(quad, u)
    t = tangents(quad, u)
    q = t[:,1] × t[:,2]
    j = norm(q)
    n = normalize(q)
    NeighborhoodQuad(quad, u, c, t, j, n)
end
function neighborhood(quad::AffineQuadrilateral, u)
    c = cartesian(quad, u)
    t = tangents(quad, u)
    NeighborhoodQuad(quad, u, c, t, quad.volume, quad.n)
end

function parametric(p::NeighborhoodQuad) p.parametric end
function cartesian(p::NeighborhoodQuad) p.cartesian end
function tangents(p::NeighborhoodQuad) p.tangents end
function tangents(p::NeighborhoodQuad, i::Int) p.tangents[:,i] end
function jacobian(p::NeighborhoodQuad) p.jacobian end
function normal(p::NeighborhoodQuad) p.normal end

function neighborhood_lazy(quad::Quadrilateral, u) neighborhood(quad, u) end


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

@testitem "Quadrilateral" begin

    #square
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)

    @test coordtype(quad) == Float64

    mp = neighborhood(quad, point(0.5, 0.5))
    @test cartesian(mp) ≈ point(0.5,0.5,0.0)
    @test tangents(mp, 1) ≈ point(1,0,0)
    @test tangents(mp, 2) ≈ point(0,1,0)
    @test jacobian(mp) ≈ 1.0

    @test quad.p21 ≈ point(1,0,0)
    @test quad.p41 ≈ point(0,1,0)
    @test quad.p1234 ≈ point(0,0,0)

    # distorted quadrilateral
    p1 = point(0.0, 0.0, 0.1)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)
    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)

    @test cartesian(neighborhood(quad, (0,0))) ≈ p1
    @test cartesian(neighborhood(quad, (1,0))) ≈ p2
    @test cartesian(neighborhood(quad, (1,1))) ≈ p3 
    @test cartesian(neighborhood(quad, (0,1))) ≈ p4 
    mp = neighborhood(quad, (0,0))
    @test tangents(mp, 1) ≈ quad.p21 ≈ p2-p1
    @test tangents(mp, 2) ≈ quad.p41 ≈ p4-p1
    @test jacobian(quad, (0.5,0.5)) > 0.0
    @test cartesian(neighborhood(quad, (0.5,0.5))) ≈ sum(quad.vertices)/4

end

@testitem "AffineQuadrilateral" begin

    # square
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)

    quad = CompScienceMeshes.AffineQuadrilateral(p1,p2,p3,p4)

    @test coordtype(quad) == Float64

    mp = neighborhood(quad, (0.5,0.5))
    @test cartesian(mp) ≈ point(0.5,0.5,0.0)
    @test tangents(mp, 1) ≈ point(1,0,0)
    @test tangents(mp, 2) ≈ point(0,1,0)
    @test quad.volume ≈ jacobian(mp) ≈ 1.0
end

@testitem "RefQuadrilateral" begin

    # distorted quadrilateral
    p1 = point(0.0, 0.0, 0.0)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)

    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)
    
    refquad = domain(quad)
    @test typeof(refquad).parameters[1] == Float64
    @test vertices(refquad)[3] == point(1, 1)

    I = [4,1,2,3]
    refquad_ = CompScienceMeshes.permute_vertices(refquad, I)

    @test refquad_.p1 ≈ point(0.0, 1.0)
    @test refquad_.a ≈ point(0.0, 0.0) - point(0.0, 1.0)
    @test refquad_.b ≈ point(1.0, 1.0) - point(0.0, 1.0) 

    @test cartesian(refquad_, (0.3,0.1)) ≈ refquad_.p1 + 0.3*refquad_.a + 0.1*refquad_.b
end



@testitem "RefQuadrilateral: quadpoints" begin
    refchart = CompScienceMeshes.RefQuadrilateral{Float64}()
    qps = quadpoints(refchart, 5)
    fn = p -> 1.0
    I = sum(w*fn(p) for (p,w) in qps)
    @test I ≈ 1.0
end

@testitem "Quadrilateral: quadpoints" begin
    p1 = point(0,0,0)
    p2 = point(2,0,0)
    p3 = point(2,3,0)
    p4 = point(0,3,0)
    quad = CompScienceMeshes.Quadrilateral(p1,p2,p3,p4)
    qps = quadpoints(quad, 5)
    fn = p -> 1.0
    I = sum(w*fn(p) for (p,w) in qps)
    @test I ≈ 6.0
end


