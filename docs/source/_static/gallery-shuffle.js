// displays a random selection of 3 gallery cards
document.addEventListener("DOMContentLoaded", () => {
    const source = document.getElementById('all-gallery-source');
    const cards = Array.from(source.querySelectorAll('.sd-col'));
    for (let i = cards.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [cards[i], cards[j]] = [cards[j], cards[i]];
    }
    cards.slice(3).forEach(card => card.remove());
    cards.slice(0, 3).forEach(card => source.querySelector('.sd-row').appendChild(card));
    source.style.display = '';
    source.removeAttribute('id');
});
