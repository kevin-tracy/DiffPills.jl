
soi = "first_term = (p2x^2 + p2y^2 + p2z^2 + 1)
second_term = (p1x^2 + p1y^2 + p1z^2 + 1)
third_term=(4*p2x^2*p2z - 8*p2x*p2y + 4*p2y^2*p2z + 4*p2z^3 - 4*p2z)
fourth_term = (4*p1x^2*p1y + 8*p1x*p1z + 4*p1y^3 + 4*p1y*p1z^2 - 4*p1y)
fifth_term = (4*p2x^2*p2y + 8*p2x*p2z + 4*p2y^3 + 4*p2y*p2z^2 - 4*p2y)
sixth_term = (8*p2y^2 + 8*p2z^2)
seventh_term = (8*p1y^2 + 8*p1z^2)
eighth_term = (p1x^4 + 2*p1x^2*p1z^2 - 8*p1x*p1y*p1z - p1y^4 + 6*p1y^2 + p1z^4 - 1)
ninth_term = (4*p1x^2*p1z - 8*p1x*p1y + 4*p1y^2*p1z + 4*p1z^3 - 4*p1z)
tenth_term = (p2x^2 - p2y^2 - p2z^2 + 1)
eleventh_term = (p2x^3 - p2x^2*p2y*p2z + p2x*p2y^2 - 3*p2x*p2z^2 + p2x - p2y^3*p2z - p2y*p2z^3 + 3*p2y*p2z)

SA[
L1*z1*(seventh_term/second_term^2 - 1) - L2*z2*(sixth_term/first_term^2 - 1)
(L1*z1*ninth_term)/second_term^2 - (L2*z2*third_term)/first_term^2
(L2*z2*fifth_term)/first_term^2 - (L1*z1*fourth_term)/second_term^2
z1*z2*((L1*L2*(8*p1y - 8*p1x*p1z)*third_term)/(second_term^2*first_term^2) - (L1*L2*(8*p1z + 8*p1x*p1y)*fifth_term)/(second_term^2*first_term^2) + (4*L1*L2*p1x*(sixth_term/first_term^2 - 1)*seventh_term)/second_term^3 + (4*L1*L2*p1x*fourth_term*fifth_term)/(second_term^3*first_term^2) + (4*L1*L2*p1x*ninth_term*third_term)/(second_term^3*first_term^2)) - z1*((L1*(8*p1z + 8*p1x*p1y)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/second_term^2 + (L1*(8*p1y - 8*p1x*p1z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/second_term^2 + (4*L1^2*fourth_term*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/second_term^5 - (4*L1^2*ninth_term*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/second_term^5 - (2*L1^2*p1x*(seventh_term/second_term^2 - 1)*seventh_term)/second_term^3 - (4*L1*p1x*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1x*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (4*L1*p1x*seventh_term*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2))/second_term^3) - z2*((4*L1*L2*third_term*(p1x^3*p1z - 3*p1x^2*p1y + p1x*p1y^2*p1z + p1x*p1z^3 - 3*p1x*p1z + p1y^3 + p1y*p1z^2 + p1y))/(second_term^3*first_term^2) - (4*L1*L2*fifth_term*(- p1x^3*p1y - 3*p1x^2*p1z - p1x*p1y^3 - p1x*p1y*p1z^2 + 3*p1x*p1y + p1y^2*p1z + p1z^3 + p1z))/(second_term^3*first_term^2) + (2*L1*L2*p1x*(sixth_term/first_term^2 - 1)*seventh_term)/second_term^3)
z2*((2*L1*L2*fifth_term*eighth_term)/(second_term^3*first_term^2) - (4*L1*L2*third_term*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/(second_term^3*first_term^2) + (8*L1*L2*p1y*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*((2*L1^2*fourth_term*eighth_term)/second_term^5 + (L1*(8*p1x - 8*p1y*p1z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/second_term^2 - (4*L1^2*ninth_term*(p1x^3 + p1x^2*p1y*p1z - 3*p1x*p1y^2 + p1x*p1z^2 + p1x + p1y^3*p1z + p1y*p1z^3 - 3*p1y*p1z))/second_term^5 + (L1*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4))/second_term^2 - (16*L1*p1y*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*p1y*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1y*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (8*L1^2*p1y*(seventh_term/second_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) + z1*z2*((L1*L2*(8*p1x - 8*p1y*p1z)*third_term)/(second_term^2*first_term^2) - (L1*L2*(4*p1x^2 + 12*p1y^2 + 4*p1z^2 - 4)*fifth_term)/(second_term^2*first_term^2) - (16*L1*L2*p1y*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 + (4*L1*L2*p1y*fourth_term*fifth_term)/(second_term^3*first_term^2) + (4*L1*L2*p1y*ninth_term*third_term)/(second_term^3*first_term^2))
z2*((2*L1*L2*third_term*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/(second_term^3*first_term^2) + (4*L1*L2*fifth_term*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/(second_term^3*first_term^2) + (8*L1*L2*p1z*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*((2*L1^2*ninth_term*(p1x^4 + 2*p1x^2*p1y^2 + 8*p1x*p1y*p1z + p1y^4 - p1z^4 + 6*p1z^2 - 1))/second_term^5 + (L1*(8*p1x + 8*p1y*p1z)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/second_term^2 + (4*L1^2*fourth_term*(p1x^3 - p1x^2*p1y*p1z + p1x*p1y^2 - 3*p1x*p1z^2 + p1x - p1y^3*p1z - p1y*p1z^3 + 3*p1y*p1z))/second_term^5 - (L1*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4))/second_term^2 - (16*L1*p1z*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*p1z*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fourth_term)/second_term^3 + (4*L1*p1z*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*ninth_term)/second_term^3 + (8*L1^2*p1z*(seventh_term/second_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3) - z1*z2*((L1*L2*(4*p1x^2 + 4*p1y^2 + 12*p1z^2 - 4)*third_term)/(second_term^2*first_term^2) + (L1*L2*(8*p1x + 8*p1y*p1z)*fifth_term)/(second_term^2*first_term^2) + (16*L1*L2*p1z*(sixth_term/first_term^2 - 1)*(p1x^2 - p1y^2 - p1z^2 + 1))/second_term^3 - (4*L1*L2*p1z*fourth_term*fifth_term)/(second_term^3*first_term^2) - (4*L1*L2*p1z*ninth_term*third_term)/(second_term^3*first_term^2))
L2*z2*(sixth_term/first_term^2 - 1) - L1*z1*(seventh_term/second_term^2 - 1)
(L2*z2*third_term)/first_term^2 - (L1*z1*ninth_term)/second_term^2
(L1*z1*fourth_term)/second_term^2 - (L2*z2*fifth_term)/first_term^2
z2*((L2*(8*p2z + 8*p2x*p2y)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/first_term^2 + (L2*(8*p2y - 8*p2x*p2z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/first_term^2 - (4*L2^2*fifth_term*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/first_term^5 + (4*L2^2*third_term*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/first_term^5 + (2*L2^2*p2x*(sixth_term/first_term^2 - 1)*sixth_term)/first_term^3 - (4*L2*p2x*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 + (4*L2*p2x*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (4*L2*p2x*sixth_term*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2))/first_term^3) - z1*((4*L1*L2*ninth_term*(p2x^3*p2z - 3*p2x^2*p2y + p2x*p2y^2*p2z + p2x*p2z^3 - 3*p2x*p2z + p2y^3 + p2y*p2z^2 + p2y))/(second_term^2*first_term^3) - (4*L1*L2*fourth_term*(- p2x^3*p2y - 3*p2x^2*p2z - p2x*p2y^3 - p2x*p2y*p2z^2 + 3*p2x*p2y + p2y^2*p2z + p2z^3 + p2z))/(second_term^2*first_term^3) + (2*L1*L2*p2x*(seventh_term/second_term^2 - 1)*sixth_term)/first_term^3) + z1*z2*((L1*L2*(8*p2y - 8*p2x*p2z)*ninth_term)/(second_term^2*first_term^2) - (L1*L2*(8*p2z + 8*p2x*p2y)*fourth_term)/(second_term^2*first_term^2) + (4*L1*L2*p2x*(seventh_term/second_term^2 - 1)*sixth_term)/first_term^3 + (4*L1*L2*p2x*fourth_term*fifth_term)/(second_term^2*first_term^3) + (4*L1*L2*p2x*ninth_term*third_term)/(second_term^2*first_term^3))
z1*((2*L1*L2*fourth_term*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/(second_term^2*first_term^3) - (4*L1*L2*ninth_term*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/(second_term^2*first_term^3) + (8*L1*L2*p2y*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3) - z2*((2*L2^2*fifth_term*(p2x^4 + 2*p2x^2*p2z^2 - 8*p2x*p2y*p2z - p2y^4 + 6*p2y^2 + p2z^4 - 1))/first_term^5 - (L2*(8*p2x - 8*p2y*p2z)*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2)))/first_term^2 - (4*L2^2*third_term*(p2x^3 + p2x^2*p2y*p2z - 3*p2x*p2y^2 + p2x*p2z^2 + p2x + p2y^3*p2z + p2y*p2z^3 - 3*p2y*p2z))/first_term^5 - (L2*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4))/first_term^2 + (16*L2*p2y*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*tenth_term)/first_term^3 + (4*L2*p2y*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 - (4*L2*p2y*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (8*L2^2*p2y*(sixth_term/first_term^2 - 1)*tenth_term)/first_term^3) + z1*z2*((L1*L2*(8*p2x - 8*p2y*p2z)*ninth_term)/(second_term^2*first_term^2) - (L1*L2*(4*p2x^2 + 12*p2y^2 + 4*p2z^2 - 4)*fourth_term)/(second_term^2*first_term^2) - (16*L1*L2*p2y*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3 + (4*L1*L2*p2y*fourth_term*fifth_term)/(second_term^2*first_term^3) + (4*L1*L2*p2y*ninth_term*third_term)/(second_term^2*first_term^3))
z1*((2*L1*L2*ninth_term*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/(second_term^2*first_term^3) + (4*L1*L2*fourth_term*eleventh_term)/(second_term^2*first_term^3) + (8*L1*L2*p2z*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3) - z2*((2*L2^2*third_term*(p2x^4 + 2*p2x^2*p2y^2 + 8*p2x*p2y*p2z + p2y^4 - p2z^4 + 6*p2z^2 - 1))/first_term^5 - (L2*(8*p2x + 8*p2y*p2z)*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2)))/first_term^2 + (4*L2^2*fifth_term*eleventh_term)/first_term^5 + (L2*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4))/first_term^2 + (16*L2*p2z*(r1x - r2x - (L1*(seventh_term/second_term^2 - 1))/2 + (L2*(sixth_term/first_term^2 - 1))/2)*tenth_term)/first_term^3 + (4*L2*p2z*(r1z - r2z + (L1*fourth_term)/(2*second_term^2) - (L2*fifth_term)/(2*first_term^2))*fifth_term)/first_term^3 - (4*L2*p2z*(r1y - r2y - (L1*ninth_term)/(2*second_term^2) + (L2*third_term)/(2*first_term^2))*third_term)/first_term^3 + (8*L2^2*p2z*(sixth_term/first_term^2 - 1)*tenth_term)/first_term^3) - z1*z2*((L1*L2*(4*p2x^2 + 4*p2y^2 + 12*p2z^2 - 4)*ninth_term)/(second_term^2*first_term^2) + (L1*L2*(8*p2x + 8*p2y*p2z)*fourth_term)/(second_term^2*first_term^2) + (16*L1*L2*p2z*(seventh_term/second_term^2 - 1)*tenth_term)/first_term^3 - (4*L1*L2*p2z*fourth_term*fifth_term)/(second_term^2*first_term^3) - (4*L1*L2*p2z*ninth_term*third_term)/(second_term^2*first_term^3))
]"

soi2 = replace(soi, "p1x" => "P1.p[0]",
                    "p1y" => "P1.p[1]",
                    "p1z" => "P1.p[2]",
                    "p2x" => "P2.p[0]",
                    "p2y" => "P2.p[1]",
                    "p2z" => "P2.p[2]",
                    "r1x" => "P1.r[0]",
                    "r1y" => "P1.r[1]",
                    "r1z" => "P1.r[2]",
                    "r2x" => "P2.r[0]",
                    "r2y" => "P2.r[1]",
                    "r2z" => "P2.r[2]",
                    "L1"  => "P1.L",
                    "L2"  => "P2.L"
                    )

# [Pair("P"*string(i)*"."*p_or_r*"["*string(idx)*"]^2", "P"*string(i)*"."*p_or_r*"["*string(idx)*"]"" * ]

L = []
for i = 1:2
    for p_or_r in ["p","r"]
        for idx = 0:2
            #Pi.p[0]
            Pi = "P"*string(i)*"."*p_or_r*"["*string(idx)*"]"
            push!(L, Pair(Pi * "^2", Pi * "*" * Pi) )
        end
    end
end
for i = 1:2
    for p_or_r in ["p","r"]
        for idx = 0:2
            #Pi.p[0]
            Pi = "P"*string(i)*"."*p_or_r*"["*string(idx)*"]"
            push!(L, Pair(Pi * "^3", "pow(" *Pi * ",3)") )
            push!(L, Pair(Pi * "^4", "pow(" *Pi * ",4)") )
            push!(L, Pair(Pi * "^5", "pow(" *Pi * ",5)") )
        end
    end
end
push!(L, Pair("P1.L^2", "P1.L*P1.L"))
push!(L, Pair("P2.L^2", "P2.L*P2.L"))
# push!(L, Pair("P1.L^3", "pow(P1.L,3)"))


for term in ["first_term","second_term","third_term","fourth_term","fifth_term","sixth_term","seventh_term","eighth_term","ninth_term","tenth_term","eleventh_term"]
    old = term * "^2"
    push!(L, Pair(old, term * "*" * term))

    old = term * "^3"
    push!(L, Pair(old, "pow(" * term * ",3)"))

    old = term * "^4"
    push!(L, Pair(old, "pow(" * term * ",4)"))

    old = term * "^5"
    push!(L, Pair(old, "pow(" * term * ",5)"))
end

soi3 = replace(soi2, L...)
