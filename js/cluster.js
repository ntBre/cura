async function addToDataset() {
	let node = document.getElementById("modal-box-content");
	let smiles = node.getAttribute("smiles");
	if (smiles) {
		let response = await fetch(`/add-molecule?smiles=${smiles}`, {
			method: "GET",
		});
		if (!response.ok) {
			console.log(`error handling request: ${result}`);
		}
	} else {
		console.log("no smiles found");
	}
}
