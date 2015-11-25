declare -i x y z

x=0
Z="{"
while [ $x -lt 18 ]
do
	Z="$Z{"
	y=0
	while [ $y -lt 31 ]
	do
		z=1
		Z="$Z{`sed -n ${z}p ${x}_${y}.pmf`"
		Z="$Z,"
		z=$z+1
		while [ $z -lt 61 ]
        	do
			Z="$Z,`sed -n ${z}p ${x}_${y}.pmf`"
			z=$z+1
        	done
		Z="$Z},"
		y=$y+1
	done
	Z="$Z},"
	x=$x+1
done
Z="$Z}"
echo $Z > test.txt
