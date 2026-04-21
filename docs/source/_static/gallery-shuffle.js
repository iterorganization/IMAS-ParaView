// displays a random selection of 3 gallery cards
document.addEventListener("DOMContentLoaded", () => {
    const source = document.getElementById('all-gallery-source');
    const cards = Array.from(source.querySelectorAll('.sd-col'));
    for (let i = 0; i < 3; i++) {
        const j = i + Math.floor(Math.random() * (cards.length - i));
        [cards[i], cards[j]] = [cards[j], cards[i]];
    }
    const row = source.querySelector('.sd-row');
    cards.slice(3).forEach(card => card.remove());
    cards.slice(0, 3).forEach(card => row.appendChild(card));
    source.style.display = '';
});
