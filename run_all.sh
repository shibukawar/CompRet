root_dir="path/to/CompRet" # Edit here
cd $root_dir
reactionLibrary=$root_dir"reactions/"
reactionBase=$root_dir
buildingBlocks=$root_dir"block.smi" # Edit here
target="OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1" # Edit here
for d in 10
do
    for t in 1800000 
    do 
        for reactions in rxns1
        do
            reactionLibrary=$reactionBase$reactions"/"
            echo $reactionLibrary
            cd $root_dir"retrosynthesis/src"
            log_dir="./log_dir"
            if [ ! -e $log_dir ]; then
                mkdir $log_dir
            fi
            log_dir=$log_dir"cetirizine_"$reactions"_"$d"_"$t"_"$routeNum"/"
            echo $log_dir
            if [ ! -e $log_dir ]; then
                mkdir $log_dir
            fi
            fname=$log_dir"cetirizine_reaxys_depth_"$d"_"$t"_"$routeNum".txt" # Edit here
            echo $fname
            if [ ! -e $fname ]; then
                echo "target:$target
depth:$d
time:$t
reaction:$reactionLibrary
building block:$buildingBlocks" > $log_dir"info.txt"
                # Modify path to ChemAxon's library
                javac -cp .:/Applications/ChemAxon/JChemSuite/lib/jchem.jar:/Applications/Junit4/junit-4.12.jar:/Applications/Junit4/hamcrest-core-2.1.jar -Xlint:unchecked,deprecation retrosynthesis/*.java
                java -cp .:/Applications/ChemAxon/JChemSuite/lib/jchem.jar:/Applications/Junit4/junit-4.12.jar:/Applications/Junit4/hamcrest-core-2.1.jar -ea retrosynthesis.Main $d $t $routeNum $reactionLibrary $buildingBlocks $target> $fname
                cp -rf $root_dir"retrosynthesis/extractDot.sh" $log_dir
                cd $log_dir
                grep "digraph" $fname > "pt.txt"
                grep "found" $fname > "timePoint.txt"
                sed -n  "/end/,/Terminate!!!/p" $fname > "finalResult.txt"
                sed -n "/StartBlocks/,/EndBlocks/p" $fname > "buildingBlockList.txt"
            fi
            cd $root_dir"extractor"
            source activate shibuchem-env
            rank=100
            tree_size=500
            tree_size=$(cat $log_dir"pt.txt" | wc -l)
            mode='enumeration'
            python main.py $reactionLibrary $log_dir $tree_size $rank $mode > $log_dir"log_"$tree_size".txt"
            route_dir=$log_dir"route"$tree_size"/"
            echo $route_dir
            echo $tree_size
            cd $route_dir
            # if DOT is installed, you can see routes in png file
            # if [ ! -e "images" ]; then
            #     mkdir "images"
            # fi
            # for k in $(ls dots)
            # do
            #     hoge=$(basename $k .dot)
            #     dot -Tpng dots/$k -o images/$hoge.png
            # done
            # cd $log_dir
            # zip -r "route"$tree_size".zip" "route"$tree_size
            #cd $root_dir"visualizer"
            #python main.py $route_dir
            #python main.py $route_dir 
        done
    done
done
