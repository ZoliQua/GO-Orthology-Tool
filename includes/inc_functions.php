<?php

/*
* ** DO NOT REMOVE **
* *******************
* Project Name: GeneOntology Extension Tool
* Project Website: http://go.orthologfindertool.com
* Project Version: Public Version 1.0
*
* Project Source Code: https://github.com/ZoliQua/GO-Extension-Tool
*
* Author: Zoltan Dul, 2018
* Email: zoltan.dul@kcl.ac.uk and zoltan.dul@gmail.com
* Twitter: @ZoliQa
*
* DESCRIPTION
* ****************
* A bioinformatics tool that aims to extend Gene Ontology and give novel suggestions for 
* funcional annotation, based on their orthological relation.
*
* PHP FILE
* *******************
* Page - MAIN File
* *******************
*
* This file redirects user handles database.
*
* *******************
*
* All code can be used under GNU General Public License version 2.
* If you have any question or find some bug please email me.
*
*/

// Classes and functions

	class VennDiagram {

		public $permutation = array();
		public $orderedPermutation = array();
		public $replacedPermutation = array();
		public $permutArray = array();
		public $translator = array();
		public $keys = array();

		public function __construct($NrOfSets = 5, $replace = false) {

			$rangeNR = range("A","M");
			$this->permuation = self::venn_permutation(range("A", $rangeNR[$NrOfSets-1] ), $NrOfSets);

			for($i = $NrOfSets;$i >= 1; $i--) $this->orderedPermutation[] = $this->permuation[$i];

			if(!$replace) $replace = array("A" => "AT", "B" => "CE", "C" => "DM", "D" => "DR", "E" => "HS", "F" => "SC", "G" => "SP");
			self::exchanger($replace);

			return $this->replacedPermutation;

		}

		private function factorial($num) {

			$returnNum = $num;

			for($i=$num-1;$i>1;$i--) $returnNum = $returnNum * $i;

			return $returnNum;
		}

		private function venn_permutation($tombom, $pos, $return = false){

			if(!$pos or $pos > count($tombom) ) return $return;

			if(!$return) $return = array();
			$return[$pos] = array();

			$level = count($tombom);

			$currentOptions = 0;

			switch ($level) {
				case $pos:
					$maxOptions = 1;
					break;
				
				default:
					$maxOptions = self::factorial($level) / self::factorial($level - $pos) / self::factorial($pos);
					break;
			}		

			while ($maxOptions != $currentOptions) {

				$rand = array_rand($tombom, $pos);
				$thisOption = array();

				if($pos != 1) foreach ($rand as $key => $v) $thisOption[] = $tombom[$v]; 
				else $thisOption[] = $tombom[$rand]; 
				
				sort($thisOption);
				$thisOption = implode(",", $thisOption);

				if(! in_array($thisOption, $return[$pos])) $return[$pos][] = $thisOption;

				$currentOptions = count($return[$pos]);			
			}

			sort($return[$pos]);
			$pos--;

			$return = self::venn_permutation($tombom, $pos, $return);

			return $return;
		}

		public function exchanger($List) {

			$newPermutation = array();
			$KeyChecker = array();
			$KeyCollector = array();
			$ConvertedList = array();

			$MockReplacR = array("A" => "XXX", "B" => "XYX", "C" => "ZYZ", "D" => "VVV", "E" => "ZXZ", "F" => "YMY", "G" => "YXY");

			foreach ($List as $from => $to) {

				$KeyCollector[$to] = array();
				$KeyChecker[] = $to;
				$ConvertedList[$MockReplacR[$from]] = $to;

			}		

			foreach ($this->orderedPermutation as $id => $stringArray) {
				
				foreach ($stringArray as $nr => $lineString) {

					$newLine = $lineString;

					foreach ($MockReplacR as $from => $to) $newLine = str_replace($from, $to, $newLine);
					foreach ($ConvertedList as $from => $to) $newLine = str_replace($from, $to, $newLine);

					$newPermutation[] = $newLine;

					foreach ($KeyChecker as $nr => $id ) {
						if(strpos($newLine, $id) !== false) $KeyCollector[$id][] = $newLine;

					}

					$this->permutArray[$newLine] = array();
					$this->translator[$newLine] = str_replace(",","", $lineString);
					
				}

			}

			$this->keys = $KeyCollector;
			$this->replacedPermutation = $newPermutation;

			return true;
		}
	}

	class SVG_File {

		public $svg = "";

		public function __construct($folders, $num, $type = "venn"){

			if($num < 2 OR $num > 7) $num = 5;

			$filename = $folders . "ortholog_" . $type . "_" . $num . "_names.svg";
			$svg_open = file_get_contents($filename);

			$this->svg = $svg_open;

			return true;
		}
	}

	function SpeciesValidation($arraySource, $species) {

		$arrayReturn = [];

		foreach ($arraySource as $k => $v) {
			if(in_array($v, $species)) $arrayReturn[] = $v;
		}

		if(count($arrayReturn) == 0) return false;
		else return $arrayReturn;
	}

	function TimeEnd($time_start, $plustxt = "Overall") {

		$time_end = microtime(true);
		$exection_time = $time_end - $time_start;

		$hours = (int) ($exection_time / 3600);
		$minutes = ( (int) ($exection_time / 60) ) - ($hours * 60);
		$seconds = $exection_time - ( ( $hours * 3600 ) + ( $minutes * 60 ) );

		$txt = $hours . " hours " . $minutes . " minutes and " . substr($seconds, 0, 5) . " seconds. [" . $exection_time . "]";

		return "<p>The <i>$plustxt</i> execution time was $txt</p>\n";
	}

	function PrintThingsOut($venn_diagram = false, $more_info = false, $num = false, $this_file, $go, $gos, $faj) {

		$mit_printeljek = "";

		if($venn_diagram) {
			$mit_printeljek .= "<div style=\"text-align: center;\">";
			$mit_printeljek .= "<H1 stlye=\"center\">".$num."-set Venn Diagram</H1>\n";
			$mit_printeljek .= "<H2 stlye=\"center\">$go - " . $gos[$go] . "</H2>";
			$mit_printeljek .= $venn_diagram;
			$mit_printeljek .= "</div>" . "<BR><BR>\n";
			$mit_printeljek .= $more_info;
		}

		$mit_printeljek .= "<BR><BR> <B>Choose a QUERY Form</B><BR>\n";
		$mit_printeljek .= "<form method='GET'>";

		$mit_printeljek .= "GO categories: <select name='thisgo'>\n";

		foreach ($gos as $k => $v) $mit_printeljek .= "<option value='$k'>$k - $v</option>\n";

		$mit_printeljek .= "</select>\n";

		$mit_printeljek .= "Species: <select name='specs[]' size=7 multiple>\n";

		foreach ($faj as $v => $arr) $mit_printeljek .= "<option value='$v'>".$arr["long"]."</option>\n";

		$mit_printeljek .= "</select>\n";
		$mit_printeljek .= " Size Manual List <INPUT type='checkbox' name='sizemanual'> ";

		$mit_printeljek .= "<INPUT type='submit' name='ok' value='Query'>\n";
		$mit_printeljek .= "</form>\n";

		return $mit_printeljek;
	}

?>