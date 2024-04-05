document.addEventListener('DOMContentLoaded', main, false);

function main() {
	let svgs = document.querySelectorAll("span");

	svgs.forEach(function (svg) {
		svg.addEventListener("click", function () {
			let smiles = svg.getAttribute("smiles");
			let natoms = svg.getAttribute("natoms");

			// get the dialog box initially in the page
			let dialog = document.getElementById("modal-box");

			// construct a new div to go into the modal, then add in the SMILES,
			// natoms, and the svg
			let frame = document.createElement("div");
			if (smiles) {
				let s = document.createElement("p");
				s.appendChild(document.createTextNode("SMILES: " + smiles));
				frame.appendChild(s);
			}
			if (natoms) {
				let n = document.createElement("p");
				n.appendChild(document.createTextNode(natoms + " atoms"));
				frame.appendChild(n);
			}
			frame.appendChild(svg.cloneNode(true));

			// put the new frame into the dialog and display it
			dialog.replaceChildren(frame);
			dialog.showModal();
		});
	});
}
