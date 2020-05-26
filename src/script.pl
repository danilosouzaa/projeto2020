use 5.010;
use POSIX qw(strftime);

sub read_lines {
    my ($file) = @_;

    open my $in, "<:encoding(utf8)", $file or die "$file: $!";
    local $/ = undef;
    my $content = <$in>;
    close $in;
    return split /\n/, $content;
}
my @instances = read_lines("../input/Instances.txt");
my @p1 =("0");# 0- Grasp, 1 - Greedy
my @p2 =("0"); # 0 - letchford, 1 -Ballas
my @p3 = ("1"); # 1 - minimal, 0- no minimal
my @p4 =("200000");## poolSize
my @p5 =("3000","4000","6000","7000", "8000","9000"); #number of iteration
my @p6 = ("-1");#alpha
my $count = 1;
my $flag = 0;
my $data = strftime "%F", localtime;
open(FIN,">>rodar.sh");
foreach my $typeSeparation (@p1){
	foreach my $typeLifted (@p2){
		if($typeLifted == 0){
			@p3 =("1");
		}else{
			@p3 = ("0","1");
		}
		foreach my $minimal (@p3){
			if($typeSeparation==1){
                        	foreach my $record (@instances) {
                                	if($record ne "\n"){
						if($flag==0){ #cria a pasta se não houver
                                       	        	print FIN ("mkdir ../outputArt/exp-$data-$count\n");
                                       	        	$flag=1;
						}
                                        	print FIN ("./projectFinal instancesGenerated/$record 10800 1000 -CC $typeSeparation $typeLifted $minimal 200000 200000 -1 >>../outputArt/exp-$data-$count/Resultado-$record.txt 2>&1 \n");
					}
				}
				$count++;
                        	$flag = 0;
			}else{
				foreach my $pool (@p4){
					foreach my $iteGrasp (@p5){
						foreach my $alpha (@p6){
                               				foreach my $record (@instances) {
                               					if($record ne "\n"){
                                        				if($flag==0){
                                        	                                print FIN ("mkdir ../outputArt/exp-$data-$count\n");
                                        	                                $flag=1;
                                        	                        }
                                        	                        print FIN ("./projectFinal instancesGenerated/$record 10800 1000 -CC $typeSeparation $typeLifted $minimal $pool $iteGrasp $alpha >>../outputArt/exp-$data-$count/Resultado-$record.txt 2>&1\n");
                                        	                }
                                        	     	}
                                        	        $count++;
                                        	        $flag = 0;
                                        	}
					}
                                }
                        }
		}
	}
}
exit;

