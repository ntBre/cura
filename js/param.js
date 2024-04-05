document.addEventListener('DOMContentLoaded', main, false);

function main() {
	let svgs = document.querySelectorAll("span");

	svgs.forEach(function (svg) {
		svg.addEventListener("click", function () {
			let smiles = svg.getAttribute("smiles");
			console.log("clicked on molecule with smiles: " + smiles);
		});
	});
}
